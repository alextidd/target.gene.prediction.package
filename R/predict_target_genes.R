#' Predict target genes of a list of fine-mapped non-coding variants for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param trait Optional. The name of the trait of interest.
#' @param tissue Optional. The tissue(s) of action for the trait.
#' @param outDir The output directory in which to save the predictions. Default is "./out".
#' @param referenceDir The directory containing the external, accompanying reference data.
#' @param variantsFile A BED file of trait-associated variants grouped by association signal, for example SNPs correlated with an index variant, or credible sets of fine-mapped variants
#' @param driversFile The file containing a list of trait driver gene symbols.
#' @param variant_to_gene_max_distance The maximum absolute distance (bp) across which variant-gene pairs are considered. Default is 2Mb. The contact data is also already filtered to 2Mb.
#'
#' @return A file of variant-gene pair predictions, with associated scores, saved in the given output directory.
#' @export
predict_target_genes <- function(trait = NULL,
                                 tissue = NULL,
                                 outDir = "out",
                                 variantsFile = "/working/lab_georgiat/alexandT/target.gene.prediction.package/external_data/reference/BC.VariantList.bed",
                                 driversFile = "/working/lab_georgiat/alexandT/target.gene.prediction.package/external_data/reference/breast_cancer_drivers_2021.txt",
                                 referenceDir = "/working/lab_georgiat/alexandT/target.gene.prediction.package/external_data/reference/",
                                 variant_to_gene_max_distance = 2e6){

  # for testing:
  # library(devtools) ; load_all() ; trait="BC" ; outDir = "out" ; variantsFile="/working/lab_georgiat/alexandT/target.gene.prediction.package/external_data/reference/BC.VariantList.bed" ; driversFile = "/working/lab_georgiat/alexandT/target.gene.prediction.package/external_data/reference/breast_cancer_drivers_2021.txt" ; referenceDir = "/working/lab_georgiat/alexandT/target.gene.prediction.package/external_data/reference/" ; variant_to_gene_max_distance = 2e6

  # silence "no visible binding" NOTE for data variables in check()
  . <- NULL

  # SETUP ======================================================================================================

  # define the output
  dir.create(outDir, recursive = T, showWarnings = F)
  out <- list()
  out$Base <- paste0(outDir,"/") %>% { if(!is.null(trait)) paste0(., trait, "_") else . }
  out$Annotations <- paste0(out$Base, "target_gene_annotations.tsv")
  out$Predictions <- paste0(out$Base, "target_gene_predictions.tsv")
  out$Performance <- paste0(out$Base, "performance.tsv")
  out$PR <- paste0(out$Base, "PrecisionRecall.pdf")

  # import the user-provided variants
  cat("Importing variants...\n")
  variants <- target.gene.prediction.package::import_BED(
    variantsFile,
    metadata_cols = c("variant", "cs"))

  # import user-provided drivers and check that all symbols are in the GENCODE data
  cat("Importing driver genes...\n")
  drivers <- target.gene.prediction.package::read_tibble(driversFile)$V1 %>%
    target.gene.prediction.package::check_driver_symbols(., driversFile)

  # import the contact data
  contact <- readRDS(paste0(referenceDir, "contact/contact.rda"))
  contact_metadata <- readRDS(paste0(referenceDir, "contact/contact_metadata.rda"))

  # import the DHSs data
  DHSs <- readRDS(paste0(referenceDir, "DHSs/DHSs.rda"))
  DHSs_metadata <- readRDS(paste0(referenceDir, "DHSs/DHSs_metadata.rda"))
  DHSs_master <-  DHSs[[1]] %>%
    dplyr::distinct(chrom, start, end, DHS) %>%
    dplyr::mutate(DHSID = paste0(".", dplyr::row_number()))
  specific_DHSs_closest_specific_genes <- readRDS(paste0(referenceDir, "DHSs/specific_DHSs_closest_specific_genes.rda"))
  specific_DHSs_closest_specific_genes_metadata <- readRDS(paste0(referenceDir, "DHSs/specific_DHSs_closest_specific_genes_metadata.rda"))

  # import the TADs data
  TADs <- readRDS(paste0(referenceDir, "TADs/TADs.rda"))

  # 1) CELL TYPE ENRICHMENT ======================================================================================================
  cat("1) Cell type enrichment...\n")

  enriched <- target.gene.prediction.package::get_enriched(DHSs$specificity,
                                                           DHSs_metadata,
                                                           contact_metadata,
                                                           specific_DHSs_closest_specific_genes_metadata,
                                                           variants)
  cat("Enriched cell type(s): ", enriched$celltypes$code, "\n")
  cat("Enriched tissue(s):", unique(enriched$celltypes$tissue), "\n")

  # Subset annotations to those in enriched tissues
  enriched$DHSs <- DHSs %>%
    lapply(dplyr::select, chrom:DHS, dplyr::contains(enriched$tissues$code))
  enriched$specific_DHSs_closest_specific_genes <- specific_DHSs_closest_specific_genes %>%
    dplyr::select(DHS, dplyr::contains(enriched$tissues$code))
  enriched$contact <- contact[dplyr::filter(enriched$tissues, !is.na(list_element))$list_element]

  # 2) ENHANCER VARIANTS ======================================================================================================
  # get variants at DHSs ('enhancer variants')
  cat("2) Finding enhancer variants...\n")
  open_variants <- DHSs_master %>%
    target.gene.prediction.package::bed_intersect_left(
      variants, .,
      keepBcoords = F,
      keepBmetadata = T)

  # The transcript-x-variant universe (masterlist of all possible transcript x open variant pairs <2Mb apart)
  txv_master <- open_variants %>%
    valr::bed_slop(both = variant_to_gene_max_distance,
                   genome = target.gene.prediction.package::ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs,
                        suffix = c(".variant", ".TSS")) %>%
    dplyr::distinct(variant = variant.variant,
                    cs = cs.variant,
                    DHS = DHS.variant,
                    enst = enst.TSS,
                    symbol = symbol.TSS
                    ) %>%
    dplyr::mutate(pair = paste0(variant, "|", enst))

  # 3) ANNOTATING ======================================================================================================
  cat("3) Annotation enhancer-gene pairs...\n")

  # 3a) ENHANCER-LEVEL INPUTS ====
  # (Intersection with list of annotations, enhancers, GWAS statistics(?), etc...)
  v <- get_v_level_annotations(open_variants,
                               DHSs,
                               variant_to_gene_max_distance)

  # 3b) TRANSCRIPT-LEVEL INPUTS ====
  t <- get_t_level_annotations(DHSs)

  # 3c) DHS-LEVEL INPUTS ====
  d <- get_d_level_annotations(open_variants)

  # 3d) CS-LEVEL INPUTS ====
  c <- get_c_level_annotations(open_variants)

  # 3e) TRANSCRIPT-X-VARIANT-LEVEL INPUTS ====
  # (3C interaction, distance, etc...)
  # The package will consider all transcripts as potential targets of a variant in CS's whose ranges are within the max distance
  txv <- get_txv_level_annotations(open_variants,
                                   variant_to_gene_max_distance,
                                   DHSs,
                                   contact,
                                   TADs)

  # 3f) GENE-X-VARIANT-LEVEL INPUTS ===
  # Summarising across transcripts within a gene
  gxv <- get_gxv_level_annotations(txv)

  # 3g) TRANSCRIPT-X-CS-LEVEL INPUTS ====
  txc <- get_txc_level_annotations(txv,
                                   open_variants)

  # 3h) TRANSCRIPT-X-DHS-LEVEL INPUTS ====
  ensts_near_vars <- unique(txv$txv_inv_distance$enst)
  txd <- get_txd_level_annotations(ensts_near_vars,
                                   DHSs_master)

  # 3i) GENE-X-DHS-LEVEL INPUTS
  gxd <- get_gxd_level_annotations(open_variants,
                                   specific_DHSs_closest_specific_genes)

  # 4) ALL INPUTS ======================================================================================================
  # Master variant-transcript matrix list
  # -> wide-format (one row per transcript|variant pair, one column per celltype)
  # -> only variant-transcript combinations within 2Mb are included
  # -> pair ID rownames: variant|enst
  cat("3) Generating master table of transcript x", trait, "variant pairs, with all annotation levels...\n")
  master <- c(v, t, d, c, txv, gxv, txc, txd, gxd) %>%
    purrr::map(~ matricise_by_pair(., txv_master))
  # master %>% write_tibble(out$Annotations)

  # MultiAssayExperiment
  sample_list <- readRDS(paste0(referenceDir, "sample_list.rda"))
  # colData
  colData <- master %>%
    lapply(colnames) %>%
    unlist %>%
    as.data.frame %>%
    dplyr::distinct() %>%
    dplyr::rename("X" = ".") %>%
    dplyr::left_join(sample_list) %>%
    tibble::column_to_rownames("X")

  # celltype-specific annotations have as many columns as there are annotated celltypes
  # non-celltype-specific annotations have one `value` column, which applies across all celltypes
  MA <- MultiAssayExperiment::MultiAssayExperiment(experiments = master, colData = colData)
  saveRDS(MA, file = paste0(out$Base, "MA.rda"))

  # 5) SCORING ======================================================================================================
  cat("5) Scoring enhancer-gene pairs...\n")
  # Generating a single score for each enhancer-gene pair, with evidence
  # sum(all non-celltype-specific values, enriched celltype-specific annotations)
  scores <- names(master) %>%
    # create comma-concatenated list of non-zero annotations contributing to the scores of each pair
    lapply(function(name)
      master[[name]] %>%
        tibble::as_tibble(rownames = "pair") %>%
        dplyr::rename_with(~ paste0(name, ".", .x), -pair) %>%
        tidyr::pivot_longer(names_to = "annotation",
                            values_to = "value",
                            cols = where(is.numeric)) %>%
        dplyr::filter(value > 0, greplany(paste0("\\.", c("value", enriched$tissues$code)), annotation))
    ) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    dplyr::group_by(pair) %>%
    dplyr::mutate(score = sum(value),
                  evidence = paste(annotation, collapse = ","),
                  n_evidence = dplyr::n_distinct(annotation)) %>%
    dplyr::ungroup() %>%
    # spread out annotations again
    tidyr::pivot_wider(id_cols = c(pair, score, evidence, n_evidence),
                       names_from = annotation,
                       values_from = value) %>%
    # get CSs and symbols
    dplyr::left_join(txv_master %>% dplyr::select(pair, cs, symbol)) %>%
    # select columns
    dplyr::select(cs, symbol, score, n_evidence, evidence, where(is.numeric))

  # 6) PERFORMANCE ======================================================================================================
  # Generate PR curves (model performance metric)
  performance <- scores %>%
    # add drivers
    dplyr::mutate(driver = symbol %in% drivers$symbol) %>%
    # get performance
    target.gene.prediction.package::get_PR(txv_master, c("score", dplyr::starts_with(names(master)))) %>%
    # add annotation level info
    purrr::map(~ dplyr::mutate(., level = sub("_.*", "", prediction_method)))

  pdf(out$PR, height = 10, onefile = T)
  performance$PR %>%
    dplyr::filter(prediction_method == "score") %>%
    target.gene.prediction.package::plot_PR(colour = prediction_type)
  performance$PR %>%
    dplyr::filter(prediction_method == "score") %>%
    dplyr::select(prediction_type, AUC) %>%
    dplyr::distinct() %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_type, AUC),
                                 y = AUC)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Predictor", y = "PR AUC")
  performance$PR %>%
    dplyr::filter(AUC == max(AUC)) %>%
    plot_PR(colour = prediction_method) +
    ggplot2::labs(x = paste("recall (n = ", unique(performance$summary$True))) +
    ggplot2::facet_wrap(~ level)
    #ggplot2::theme(legend.position = "none")
  performance$PR %>%
    dplyr::select(prediction_method, prediction_type, AUC) %>%
    dplyr::filter(prediction_type == "score") %>%
    dplyr::distinct() %>%
    dplyr::mutate(level = prediction_method %>% gsub("_.*", "", .)) %>%
    dplyr::distinct() %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_method, AUC),
                                 y = AUC,
                                 fill = level)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Predictor", y = "PR AUC") +
    ggplot2::coord_flip()
  dev.off()

  # 7) XGBoost MODEL TRAINING ======================================================================================================
  # format training set
  full <- names(master) %>%
    lapply(function(name)
      master[[name]] %>%
        tibble::as_tibble(rownames = "pair") %>%
        dplyr::rename_with(~ paste0(name, ".", .x), -pair)
      ) %>%
    purrr::reduce(dplyr::full_join, by = "pair") %>%
    dplyr::mutate(label = greplany(drivers$enst, pair) %>% as.numeric)
  train <- list(data = full %>% dplyr::select(-c(pair, label)) %>% as.matrix,
                label = full$label)
  dtrain <- xgboost::xgb.DMatrix(data = train$data ,
                                 label = train$label)
  # model training
  xgb1 <- xgboost::xgboost(data = dtrain,
                           max.depth = 2,
                           eta = 1,
                           nrounds = 100,
                           objective = "binary:logistic",
                           verbose = 1)
  # view feature importance plot
  xgb1_feature_importance_mat <- xgboost::xgb.importance(feature_names = colnames(dtrain), model = xgb1)
  xgb1_feature_importance_plot <- xgboost::xgb.ggplot.importance(importance_matrix = xgb1_feature_importance_mat, top_n=50)
  # save plot
  pdf(paste0(out$Base, "xgb_model_feature_importance.pdf"))
  print(xgb1_feature_importance_plot)
  dev.off()

  return(MA)
}




# # scores for final cs|gene predictions
# gxc_predictions <- txv_scores %>%
#   dplyr::left_join(txv_master %>% dplyr::select(pair, symbol, cs)) %>%
#   # filter to the maximum-scoring pair for each CS-gene combination
#   dplyr::group_by(cs, symbol) %>%
#   dplyr::filter(score == max(score)) %>%
#   # find the maximum-scoring gene for each CS
#   dplyr::group_by(cs) %>%
#   dplyr::mutate(max_score = as.numeric(score == max(score))) %>%
#   # find predictions (score > mean(score))
#   dplyr::ungroup() %>%
#   dplyr::mutate(prediction = as.numeric(score > mean(score))) %>%
#   # calculate number of pieces of pair-supporting evidence
#   dplyr::mutate(n_evidence = stringr::str_count(evidence, ",") + 1) %>%
#   # finalise columns
#   dplyr::select(cs, symbol, score, prediction, max_score, n_evidence, evidence)
# # Generate PR curves (model performance metric)
# performance <- target.gene.prediction.package::get_PR(predictions, score)
# # PR of each individual annotation (columns of master)
# performance_all <- target.gene.prediction.package::get_PR(master, dplyr::all_of(annotation_cols))
# # add annotation level info
# performance_all$PR <- performance_all$PR %>% dplyr::mutate(level = sub("_.*", "", prediction_type))
#
# pdf(out$PR, height = 10, onefile = T)
# performance$PR %>% target.gene.prediction.package::plot_PR(colour = prediction_type)
# performance$PR %>%
#   dplyr::select(prediction_type, AUC) %>%
#   dplyr::distinct() %>%
#   ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_type, AUC),
#                                y = AUC)) +
#   ggplot2::geom_col() +
#   ggplot2::coord_flip() +
#   ggplot2::labs(x = "Predictor", y = "PR AUC")
# performance_all$PR %>%
#   dplyr::filter(AUC == max(AUC)) %>%
#   plot_PR(colour = prediction_method) +
#   ggplot2::theme(legend.position = "none")
# performance_all$PR %>%
#   dplyr::select(prediction_method, prediction_type, AUC) %>%
#   dplyr::filter(prediction_type == "score") %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(level = prediction_method %>% gsub("_.*", "", .)) %>%
#   dplyr::distinct() %>%
#   ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_method, AUC),
#                                y = AUC,
#                                fill = level)) +
#   ggplot2::geom_col() +
#   ggplot2::labs(x = "Predictor", y = "PR AUC") +
#   ggplot2::facet_wrap( ~ level, scales = "free_x")
# dev.off()


