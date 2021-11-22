#' Predict target genes of fine-mapped variants for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param trait Optional. The name of the trait of interest.
#' @param tissue_of_interest Optional. The tissue(s) of interest for the trait.
#' @param outDir The output directory in which to save the predictions. Default is "./out".
#' @param variantsFile A BED file of trait-associated variants grouped by association signal, for example SNPs correlated with an index variant, or credible sets of fine-mapped variants
#' @param driversFile The file containing a list of trait driver gene symbols.
#' @param referenceDir The directory containing the external, accompanying reference data.
#' @param variant_to_gene_max_distance The maximum absolute distance (bp) across which variant-gene pairs are considered. Default is 2Mb. The contact data is also already filtered to 2Mb.
#' @param min_proportion_of_variants_in_top_DHSs A threshold proportion of variants that reside in the specific DHSs of a celltype for that celltype to be considered enriched. Default is 5% (0.05).
#' @param include_all_celltypes_in_the_enriched_tissue If TRUE, the package will combine annotations across all available cell types within the enriched tissue, not just those in the exact enriched cell type.
#' @param do_all_cells If TRUE, the package will combine annotations across all available cell types, not just those with enriched enhancer variants. Default is FALSE.
#' @param do_manual_weighting If TRUE, runs the manual weighting chunk of the script (weight_and_score_manually) to test out different combinations of annotations to generate a score. Default is FALSE.
#' @param n_unique_manual_weights The number of unique weights for the do_manual_weighting chunk to consider. If NULL, the chunk will consider as many unique weights as there are to_add components. Default is NULL.
#' @param do_scoring If TRUE, runs the scoring chunk of the script, which combines all of the constituent MAE annotations into one score per transcript-variant pair. Default is FALSE.
#' @param do_performance If TRUE, runs the performance chunk of the script, which measures the performance of the score and each of its constituent annotations in predicting drivers as the targets of nearby variants. Default is FALSE.
#' @param do_XGBoost If TRUE, runs the XGBoost chunk of the script, which generates a model to predict the targets of variants from all available annotations and rates the importance of each annotation. Default is FALSE.
#' @param contact If you are repeatedly running predict_target_genes, you can load the `contact` object from the referenceDir into the global environment and pass it to the function to prevent redundant re-loading each call to predict_target_genes.
#' @param DHSs If you are repeatedly running predict_target_genes, you can load the `DHSs` object from the referenceDir into the global environment and pass it to the function to prevent redundant re-loading with each call to predict_target_genes.
#' @return A MultiAssayExperiment object with one assay object per annotation, one row per variant-transcript pair and one column per cell type (or 'value' if it is a non-cell-type-specific annotation).
#' @export
predict_target_genes <- function(trait = NULL,
                                 tissue_of_interest = NULL,
                                 outDir = "out",
                                 variantsFile = "/working/lab_georgiat/alexandT/tgp/example_data/data/BC.VariantList.bed",
                                 driversFile = "/working/lab_georgiat/alexandT/tgp/example_data/data/breast_cancer_drivers_2021.txt",
                                 referenceDir = "/working/lab_georgiat/alexandT/tgp/reference_data/data/",
                                 variant_to_gene_max_distance = 2e6,
                                 min_proportion_of_variants_in_top_DHSs = 0.05,
                                 include_all_celltypes_in_the_enriched_tissue = T,
                                 do_all_cells = F,
                                 do_manual_weighting = F,
                                 n_unique_manual_weights = NULL,
                                 do_scoring = F,
                                 do_performance = F,
                                 do_XGBoost = F,
                                 contact = NULL,
                                 DHSs = NULL){

  # for testing internally:
  # setwd("/working/lab_georgiat/alexandT/tgp") ; library(devtools) ; load_all() ; tissue_of_interest = NULL ; trait="BC" ; outDir = "out/BC_enriched_cells/" ; variantsFile="/working/lab_georgiat/alexandT/tgp/example_data/data/BC.VariantList.bed" ; driversFile = "/working/lab_georgiat/alexandT/tgp/example_data/data/breast_cancer_drivers_2021.txt" ; referenceDir = "/working/lab_georgiat/alexandT/tgp/reference_data/data/" ; variant_to_gene_max_distance = 2e6 ; min_proportion_of_variants_in_top_DHSs = 0.05 ; include_all_celltypes_in_the_enriched_tissue = T ; do_all_cells = F ; do_manual_weighting = F ; n_unique_manual_weights = NULL ; do_scoring = T ; do_performance = T ; do_XGBoost = T ;
  # contact = NULL ; DHSs = NULL
  # for testing externally:
  # library(devtools) ; setwd("/working/lab_georgiat/alexandT/tgp") ; load_all() ; referenceDir = "/working/lab_georgiat/alexandT/tgp/reference_data/data/" ; DHSs <- readRDS(paste0(referenceDir, "DHSs/DHSs.rda")) ; contact <- readRDS(paste0(referenceDir, "contact/contact.rda")) ; MA <- predict_target_genes(outDir = "out/BC_enriched_cells/", include_all_celltypes_in_the_enriched_tissue = F, contact = contact, DHSs = DHSs)

  # silence "no visible binding" NOTE for data variables in check()
  . <- NULL

  # SETUP ======================================================================================================
  if(do_XGBoost){do_scoring <- T}

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
  variants <- import_BED(
    variantsFile,
    metadata_cols = c("variant", "cs"))

  # import user-provided drivers and check that all symbols are in the GENCODE data
  cat("Importing driver genes...\n")
  drivers <- read_tibble(driversFile)$V1 %>%
    check_driver_symbols(., driversFile)

  # import the contact data
  cat("Importing contact data...\n")
  if(is.null(contact)){contact <- readRDS(paste0(referenceDir, "contact/contact.rda"))}

  # import the DHS binning data
  cat("Importing DHS binning data...\n")
  if(is.null(DHSs)){DHSs <- readRDS(paste0(referenceDir, "DHSs/DHSs.rda"))}
  DHSs_master <- DHSs[[1]] %>%
    dplyr::distinct(chrom, start, end, DHS)
  specific_DHSs_closest_specific_genes <- readRDS(paste0(referenceDir, "DHSs/specific_DHSs_closest_specific_genes.rda"))

  # import the expression data
  cat("Importing RNA-seq expression data...\n")
  expression <- read.delim(paste0(referenceDir, "expression.tsv"))

  # import the TADs data
  cat("Importing TAD data...\n")
  TADs <- readRDS(paste0(referenceDir, "TADs/TADs.rda"))

  # metadata for all annotations
  all_metadata <- read_tibble(paste0(referenceDir, "all_metadata.tsv"), header = T)

  # 1) CELL TYPE ENRICHMENT ======================================================================================================
  cat("1) Cell type enrichment...\n")
  enriched <- get_enriched(variants,
                           DHSs,
                           specific_DHSs_closest_specific_genes,
                           contact,
                           expression,
                           TADs,
                           all_metadata,
                           min_proportion_of_variants_in_top_DHSs,
                           tissue_of_interest,
                           do_all_cells,
                           include_all_celltypes_in_the_enriched_tissue)

  # 2) ENHANCER VARIANTS ======================================================================================================
  # get variants at DHSs ('enhancer variants')
  cat("2) Finding enhancer variants...\n")
  open_variants <- DHSs_master %>%
    bed_intersect_left(
      variants, .,
      keepBcoords = F,
      keepBmetadata = F)

  # The transcript-x-variant universe (masterlist of all possible transcript x variant pairs <2Mb apart)
  txv_master <- variants %>%
    valr::bed_slop(both = variant_to_gene_max_distance,
                   genome = ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., TSSs,
                        suffix = c(".variant", ".TSS")) %>%
    dplyr::distinct(chrom,
                    start.variant = start.variant + variant_to_gene_max_distance,
                    end.variant = end.variant - variant_to_gene_max_distance,
                    start.TSS, end.TSS,
                    variant = variant.variant,
                    pair = paste0(variant.variant, "|", enst.TSS),
                    cs = cs.variant,
                    enst = enst.TSS,
                    symbol = symbol.TSS,
                    ensg = ensg.TSS)

  # 3) ANNOTATING ======================================================================================================
  cat("3) Annotating enhancer-gene pairs...\n")

  # 3a) VARIANT-LEVEL INPUTS ====
  v <- get_v_level_annotations(variants,
                               enriched,
                               txv_master)

  # 3b) TRANSCRIPT-LEVEL INPUTS ====
  t <- get_t_level_annotations(TSSs,
                               enriched)

  # 3c) GENE-LEVEL INPUTS ===
  g <- get_g_level_annotations(txv_master,
                               enriched)

  # 3c) CS-LEVEL INPUTS ====
  c <- get_c_level_annotations(variants)

  # 3d) TRANSCRIPT-X-VARIANT-LEVEL INPUTS ====
  txv <- get_txv_level_annotations(variants,
                                   txv_master,
                                   variant_to_gene_max_distance,
                                   enriched)

  # 3e) GENE-X-VARIANT-LEVEL INPUTS ===
  gxv <- get_gxv_level_annotations(txv,
                                   variants,
                                   DHSs_master,
                                   enriched)

  # 3f) TRANSCRIPT-X-CS-LEVEL INPUTS ====
  txc <- get_txc_level_annotations(txv,
                                   variants)

  # 4) ALL INPUTS ======================================================================================================
  # Master variant-transcript matrix list
  # -> wide-format (one row per transcript|variant pair, one column per celltype)
  # -> only variant-transcript combinations within 2Mb are included
  # -> pair ID rownames: variant|enst
  cat("3) Generating master table of transcript x", trait, "variant pairs, with all annotation levels...\n")
  master <- c(v, t, g, c, txv, gxv, txc) %>%
    purrr::map(~ matricise_by_pair(., txv_master))
  # master %>% write_tibble(out$Annotations)

  # MultiAssayExperiment
  # colData
  colData <- master %>%
    lapply(colnames) %>%
    unlist %>%
    as.data.frame %>%
    dplyr::distinct() %>%
    dplyr::rename(name = ".") %>%
    dplyr::left_join(all_metadata %>% dplyr::distinct(name, tissue)) %>%
    tibble::column_to_rownames("name")

  # celltype-specific annotations have as many columns as there are annotated celltypes
  # non-celltype-specific annotations have one `value` column, which applies across all celltypes
  MA <- MultiAssayExperiment::MultiAssayExperiment(experiments = master, colData = colData)
  saveRDS(MA, file = paste0(out$Base, "MA.rda"))

  # 5) WEIGHTING ======================================================================================================
  # MANUAL WEIGHTING ===
  if(do_manual_weighting){
  celltype_of_interest <- unique(enriched$celltypes$name)
  if(length(celltype_of_interest) == 1){
    to_add <- c("v_DHSs_signal",
                "v_DHSs_specificity",
                "txv_contact_ChIAPET_binary",
                "txc_n_multicontact_binary_ChIAPET",
                "gxv_specific_DHSs_closest_specific_genes",
                "txv_exon")
    to_multiply <- c("txv_TADs",
                     "g_expression")
    manual_models <- weight_and_score_manually(MA,
                                               celltype_of_interest,
                                               txv_master,
                                               drivers,
                                               to_add,
                                               to_multiply,
                                               n_unique_manual_weights)
    write_tibble(manual_models, paste0(out$Base, "manual_weighting_models_performance.tsv"))
  } else { cat("dplyr::n_distinct(enriched$celltypes$name) != 1\nFunction will not work\n") }
  }

  # 6) SCORING ======================================================================================================
  if(do_scoring){
  cat("5) Scoring enhancer-gene pairs...\n")
  # Generating a single score for each variant-transcript pair, with evidence
  # sum(all non-celltype-specific values, enriched celltype-specific annotations)
  scores <- names(master) %>%
    # create comma-concatenated list of non-zero annotations contributing to the scores of each pair
    lapply(function(name)
      master[[name]] %>%
        tibble::as_tibble(rownames = "pair") %>%
        dplyr::rename_with(~ paste0(name, ".", .x), -pair) %>%
        tidyr::pivot_longer(names_to = "annotation",
                            values_to = "value",
                            cols = where(is.numeric)) #%>%
        #dplyr::filter(value > 0, greplany(paste0("\\.", c("value", enriched$tissues$code)), annotation))
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

  # predictions to save (max gene per CS)
  predictions <- scores %>%
    # filter to max score per CS-gene pair
    dplyr::group_by(cs, symbol) %>%
    dplyr::filter(score == max(score)) %>%
    # annotate max gene per CS
    dplyr::group_by(cs) %>%
    dplyr::mutate(max = score == max(score)) %>%
    # order
    dplyr::select(cs, symbol, score, max) %>%
    dplyr::arrange(-score)

  # write tables
  write_tibble(scores, filename = out$Annotations)
  write_tibble(predictions, filename = out$Predictions)
  }

  # 7) PERFORMANCE ======================================================================================================
  if(do_performance){
  # Generate PR curves (model performance metric)
  performance <- scores %>%
    # get performance
    get_PR(txv_master, drivers, c("score", dplyr::starts_with(names(master)))) %>%
    # add annotation level info
    purrr::map(~ dplyr::mutate(., level = sub("_.*", "", prediction_method)))

  pdf(out$PR, height = 10, width = 20, onefile = T)
  performance$PR %>%
    dplyr::filter(prediction_method == "score") %>%
    plot_PR(colour = prediction_type)
  performance$PR %>%
    dplyr::filter(prediction_method == "score") %>%
    dplyr::select(prediction_type, PR_AUC) %>%
    dplyr::distinct() %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_type, PR_AUC),
                                 y = PR_AUC)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Predictor", y = "PR AUC")
  performance$PR %>%
    dplyr::filter(PR_AUC == max(PR_AUC)) %>%
    plot_PR(colour = prediction_method) +
    # ggplot2::geom_label(data = . %>%
    #                      dplyr::ungroup() %>%
    #                      dplyr::filter(PR_AUC == max(PR_AUC),
    #                                    prediction_type == "max",
    #                                    .threshold %ni% c(Inf, 0)),
    #                    ggplot2::aes(label = prediction_method)) +
    ggplot2::labs(x = paste0("recall (n = ", unique(performance$summary$True), ")"))
  performance$PR %>%
    dplyr::select(prediction_method, prediction_type, PR_AUC) %>%
    dplyr::filter(prediction_type == "score") %>%
    dplyr::distinct() %>%
    dplyr::mutate(level = prediction_method %>% gsub("_.*", "", .)) %>%
    dplyr::distinct() %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_method, PR_AUC),
                                 y = PR_AUC,
                                 fill = level)) +
    ggplot2::geom_col() +
    ggplot2::labs(x = "Predictor", y = "PR AUC") +
    ggplot2::coord_flip()
  performance$summary %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_method, Fscore),
                                 y = Fscore,
                                 fill = level)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(~prediction_type) +
    ggplot2::coord_flip()
  dev.off()
  }

  # 8) XGBoost MODEL TRAINING ======================================================================================================
  if(do_XGBoost){
  # format training set
  full <- scores %>%
    dplyr::mutate(label = (symbol %in% drivers$symbol) %>% as.numeric) %>%
    dplyr::group_by(cs) %>%
    dplyr::filter(any(label == T)) %>%
    dplyr::ungroup()
  train <- list(data = full %>% dplyr::select(-c(cs, symbol, evidence, n_evidence, label)) %>% as.matrix,
                label = full$label)
  dtrain <- xgboost::xgb.DMatrix(data = train$data,
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
  }

  # 9) SAVE ===
  # save(master,
  #      predictions,
  #      performance,
  #      xgb1,
  #      file = paste0(out$Base, ".Rdata"))

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
# performance <- get_PR(predictions, score)
# # PR of each individual annotation (columns of master)
# performance_all <- get_PR(master, dplyr::all_of(annotation_cols))
# # add annotation level info
# performance_all$PR <- performance_all$PR %>% dplyr::mutate(level = sub("_.*", "", prediction_type))
#
# pdf(out$PR, height = 10, onefile = T)
# performance$PR %>% plot_PR(colour = prediction_type)
# performance$PR %>%
#   dplyr::select(prediction_type, PR_AUC) %>%
#   dplyr::distinct() %>%
#   ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_type, PR_AUC),
#                                y = PR_AUC)) +
#   ggplot2::geom_col() +
#   ggplot2::coord_flip() +
#   ggplot2::labs(x = "Predictor", y = "PR PR_AUC")
# performance_all$PR %>%
#   dplyr::filter(PR_AUC == max(PR_AUC)) %>%
#   plot_PR(colour = prediction_method) +
#   ggplot2::theme(legend.position = "none")
# performance_all$PR %>%
#   dplyr::select(prediction_method, prediction_type, PR_AUC) %>%
#   dplyr::filter(prediction_type == "score") %>%
#   dplyr::distinct() %>%
#   dplyr::mutate(level = prediction_method %>% gsub("_.*", "", .)) %>%
#   dplyr::distinct() %>%
#   ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_method, PR_AUC),
#                                y = PR_AUC,
#                                fill = level)) +
#   ggplot2::geom_col() +
#   ggplot2::labs(x = "Predictor", y = "PR PR_AUC") +
#   ggplot2::facet_wrap( ~ level, scales = "free_x")
# dev.off()


