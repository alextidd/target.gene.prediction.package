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
#' @param do_all_cells If TRUE, the package will combine annotations across all available cell types, not just those with enriched enhancer variants. Default is FALSE.
#' @param do_all_cells_in_enriched_tissue If TRUE, the package will combine all annotations for the tissue of the enriched celltype(s), not just the specifically enriched celltype(s). Default is TRUE. Make sure the enriched celltype has coverage across all annotations (TADs, contact, expression, DHSs) in the metadata table.
#' @param do_manual_weighting If TRUE, runs the manual weighting chunk of the script (weight_and_score_manually) to test out different combinations of annotations to generate a score. Default is FALSE.
#' @param n_unique_manual_weights The number of unique weights for the do_manual_weighting chunk to consider. If NULL, the chunk will consider as many unique weights as there are to_add components. Default is NULL.
#' @param do_scoring If TRUE, runs the scoring chunk of the script, which combines all of the constituent MAE annotations into one score per transcript-variant pair. Default is FALSE.
#' @param do_performance If TRUE, runs the performance chunk of the script, which measures the performance of the score and each of its constituent annotations in predicting drivers as the targets of nearby variants. Default is FALSE.
#' @param do_XGBoost If TRUE, runs the XGBoost chunk of the script, which generates a model to predict the targets of variants from all available annotations and rates the importance of each annotation. Default is FALSE.
#' @param do_timestamp If TRUE, will save output into a subdirectory timestamped with the data/time of the run.
#' @param contact If you are repeatedly running predict_target_genes, you can load the `contact` object from the referenceDir into the global environment and pass it to the function to prevent redundant re-loading each call to predict_target_genes.
#' @param DHSs If you are repeatedly running predict_target_genes, you can load the `DHSs` object from the referenceDir into the global environment and pass it to the function to prevent redundant re-loading with each call to predict_target_genes.
#' @return A MultiAssayExperiment object with one assay object per annotation, one row per variant-transcript pair and one column per cell type (or 'value' if it is a non-cell-type-specific annotation).
#' @export
predict_target_genes <- function(trait = NULL,
                                 tissue_of_interest = NULL,
                                 outDir = "out",
                                 variantsFile = "/working/lab_jonathb/alexandT/tgp/example_data/data/BC/BC.VariantList.bed",
                                 driversFile = "/working/lab_jonathb/alexandT/tgp/example_data/data/BC/breast_cancer_drivers_2021.txt",
                                 referenceDir = "/working/lab_jonathb/alexandT/tgp/reference_data/data/",
                                 variant_to_gene_max_distance = 2e6,
                                 min_proportion_of_variants_in_top_DHSs = 0.05,
                                 do_all_cells = F,
                                 do_all_cells_in_enriched_tissue = T,
                                 do_manual_weighting = F,
                                 n_unique_manual_weights = NULL,
                                 do_scoring = T,
                                 do_performance = T,
                                 do_XGBoost = F,
                                 do_timestamp = F,
                                 contact = NULL,
                                 DHSs = NULL){

  # for testing internally:
  # contact = NULL ; DHSs = NULL ; setwd("/working/lab_jonathb/alexandT/tgp") ; library(devtools) ; load_all() ; tissue_of_interest = NULL ; trait="BC" ; outDir = "out/" ; variantsFile="/working/lab_jonathb/alexandT/tgp/example_data/data/BC/BC.VariantList.bed" ; driversFile = "/working/lab_jonathb/alexandT/tgp/example_data/data/BC/BC.Drivers.txt" ; referenceDir = "/working/lab_jonathb/alexandT/tgp/reference_data/data/" ; variant_to_gene_max_distance = 2e6 ; min_proportion_of_variants_in_top_DHSs = 0.05 ; do_all_cells = F ; do_all_cells_in_enriched_tissue = T ; do_manual_weighting = F ; n_unique_manual_weights = NULL ; do_scoring = T ; do_performance = T ; do_XGBoost = T ; do_timestamp = F ;
  # for testing externally:
  # library(devtools) ; setwd("/working/lab_jonathb/alexandT/tgp") ; load_all() ; referenceDir = "/working/lab_jonathb/alexandT/tgp/reference_data/data/" ; DHSs <- readRDS(paste0(referenceDir, "DHSs/DHSs.rda")) ; contact <- readRDS(paste0(referenceDir, "contact/contact.rda")) ; MA <- predict_target_genes(outDir = "out/BC_enriched_cells/", contact = contact, DHSs = DHSs)

  # silence "no visible binding" NOTE for data variables in check()
  . <- NULL

  # SETUP ======================================================================================================
  if(do_XGBoost){do_scoring <- T}

  # define the output
  out <- list(
    Base = "",
    TissueEnrichment = "tissue_fisher_enrichments.tsv",
    Annotations = "target_gene_annotations.tsv",
    Predictions = "target_gene_predictions_full.tsv",
    MaxPredictions = "target_gene_predictions_max.tsv",
    Performance = "performance.tsv",
    PR = "PrecisionRecall.pdf"
  ) ; out <- paste0(outDir,"/", trait, "/") %>%
  { if(do_all_cells) paste0(., "all_") else paste0(., "enriched_") } %>%
  { if(do_all_cells_in_enriched_tissue) paste0(., "tissues/") else paste0(., "cells/") } %>%
  { if(do_timestamp) paste0(., format(Sys.time(), "%Y%m%d_%H%M%S"), "/") else . } %>%
  { purrr::map(out, function(x) paste0(., x)) }
  dir.create(out$Base, showWarnings = F, recursive = T)

  # import the user-provided variants
  cat(" > Importing variants...\n")
  variants <- import_BED(
    variantsFile,
    metadata_cols = c("variant", "cs"))

  # import user-provided drivers and check that all symbols are in the GENCODE database
  cat(" > Importing driver genes...\n")
  drivers <- read_tibble(driversFile)$V1 %>%
    check_driver_symbols(., driversFile)

  # import the contact data
  if(is.null(contact)){
    cat(" > Importing contact data...\n")
    contact <- readRDS(paste0(referenceDir, "contact.rda"))
  }

  # import the DHS binning data
  if(is.null(DHSs)){
    cat(" > Importing DHS binning data...\n")
    DHSs <- readRDS(paste0(referenceDir, "DHSs.rda"))
  }
  specific_DHSs_closest_specific_genes <- readRDS(paste0(referenceDir, "specific_DHSs_closest_specific_genes.rda"))

  # import the expression data
  cat(" > Importing RNA-seq expression data...\n")
  expression <- read.delim(paste0(referenceDir, "expression.tsv"))

  # import the TADs data
  cat(" > Importing TAD data...\n")
  TADs <- readRDS(paste0(referenceDir, "TADs.rda"))

  # metadata for all annotations
  all_metadata <- read_tibble(paste0(referenceDir, "all_metadata.tsv"), header = T)

  # 1) CELL TYPE ENRICHMENT ======================================================================================================
  cat("1) Performing cell type enrichment...\n")
  enriched <- get_enriched(variants,
                           DHSs,
                           specific_DHSs_closest_specific_genes,
                           contact,
                           expression,
                           TADs,
                           all_metadata,
                           out,
                           min_proportion_of_variants_in_top_DHSs,
                           tissue_of_interest,
                           do_all_cells,
                           do_all_cells_in_enriched_tissue)

  # 2) ENHANCER VARIANTS ======================================================================================================
  # get variants at DHSs ('enhancer variants')
  cat("2) Finding enhancer variants...\n")

  # The transcript-x-variant universe (masterlist of all possible transcript x variant pairs < variant_to_gene_max_distance apart)
  nearby_genes <- variants %>%
    # all genes within variant_to_gene_max_distance
    valr::bed_slop(both = variant_to_gene_max_distance,
                   genome = ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., TSSs,
                        suffix = c("", ".TSS")) %>%
    # fix coords
    dplyr::select(-c(start, -end))

  txv_master <- variants %>%
    # all genes within 2Mb
    valr::bed_slop(both = variant_to_gene_max_distance,
                   genome = ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., TSSs,
                        suffix = c("", ".TSS")) %>%
    # restore coords
    dplyr::select(-c(start, end)) %>%
    dplyr::left_join(variants, by = c("chrom", "variant", "cs")) %>%
    dplyr::transmute(chrom,
                    start.variant = start,
                    end.variant = end,
                    start.TSS, end.TSS,
                    distance = abs(end - end.TSS),
                    variant,
                    pair = paste0(variant, "|", enst.TSS),
                    cs,
                    enst = enst.TSS,
                    symbol = symbol.TSS,
                    ensg = ensg.TSS) %>%
    dplyr::distinct()

  # 3) ANNOTATING ======================================================================================================
  cat("3) Annotating enhancer-gene pairs...\n")

  # 3a) VARIANT-LEVEL INPUTS ====
  cat(" > V\tAnnotating variants...\n")
  v <- get_v_level_annotations(variants,
                               DHSs,
                               enriched,
                               txv_master)

  # 3b) TRANSCRIPT-LEVEL INPUTS ====
  cat(" > T\tAnnotating transcripts...\n")
  t <- get_t_level_annotations(TSSs,
                               enriched)

  # 3c) GENE-LEVEL INPUTS ===
  cat(" > G\tAnnotating genes...\n")
  g <- get_g_level_annotations(txv_master,
                               enriched)

  # 3d) CS-LEVEL INPUTS ====
  cat(" > C\tAnnotating credible sets...\n")
  c <- get_c_level_annotations(variants)

  # 3e) TRANSCRIPT-X-VARIANT-LEVEL INPUTS ====
  cat(" > TxV\tAnnotating transcript x variant pairs...\n")
  txv <- get_txv_level_annotations(variants,
                                   txv_master,
                                   variant_to_gene_max_distance,
                                   enriched)

  # 3f) GENE-X-VARIANT-LEVEL INPUTS ===
  cat(" > GxV\tAnnotating gene x variant pairs...\n")
  gxv <- get_gxv_level_annotations(variants,
                                   txv_master,
                                   DHSs_master,
                                   enriched)

  # 3g) TRANSCRIPT-X-CS-LEVEL INPUTS ====
  cat(" > TxC\tAnnotating transcript x credible set pairs...\n")
  txc <- get_txc_level_annotations(txv,
                                   variants)

  # 4) ALL INPUTS ======================================================================================================
  # Master variant-transcript matrix list
  # -> wide-format (one row per transcript|variant pair, one column per celltype)
  # -> only variant-transcript combinations within 2Mb are included
  # -> pair ID rownames: variant|enst
  # Each annotation will be aggregated across samples, taking the maximum value per pair
  cat("4) Generating master table of transcript x", trait, "variant pairs, with all annotation levels...\n")
  master <- c(v, t, g, c, txv, gxv, txc) %>%
    purrr::map(~ matricise_by_pair(., txv_master))

  # MultiAssayExperiment colData
  colData <- master %>%
    lapply(colnames) %>%
    unlist %>%
    as.data.frame %>%
    dplyr::distinct() %>%
    dplyr::rename(name = ".") %>%
    dplyr::left_join(all_metadata %>% dplyr::distinct(name, tissue), by = "name") %>%
    tibble::column_to_rownames("name")

  # celltype-specific annotations have as many columns as there are annotated celltypes
  # non-celltype-specific annotations have one `value` column, which applies across all celltypes
  MA <- MultiAssayExperiment::MultiAssayExperiment(experiments = master, colData = colData)
  saveRDS(MA, file = paste0(out$Base, "MA.rda"))

  # 5) SCORING ======================================================================================================
  if(do_scoring){
  cat("5) Scoring enhancer-gene pairs...\n")

  # declare weights (ALL OTHERS = 0.33)
  weights <- list(
    txv_TADs = 1,
    txv_inv_distance = 1,
    gxv_specific_DHSs_closest_specific_genes = 1,
    txv_intron = 1,
    gxv_missense = 1,
    gxv_nonsense = 1,
    gxv_splicesite = 1,
    t_DHSs_signal = 1,
    g_expression = 1,
    v_inv_n_genes = 0.66,
    c_inv_n_variants = 0.66
  )

  # Generating a single score for each variant-transcript pair, with evidence
  # mean(all non-celltype-specific values, enriched celltype-specific annotations) * (expression binary)
  scores <- master %>% names %>%
    sapply(function(name) {
      # weight annotations
      master[[name]] %>%
        tibble::as_tibble(rownames = "pair") %>%
        dplyr::mutate(
          value_unweighted = value,
          value = value * { if (name %in% names(weights)) weights[[name]] else 0.33 }
        )
    }, USE.NAMES = T, simplify = F) %>%
    dplyr::bind_rows(.id = "prediction_method") %>%
    # calculate score (mean annotation value)
    dplyr::group_by(pair) %>%
    dplyr::mutate(
      score = mean(value),
      score_unweighted = mean(value_unweighted)
    ) %>%
    dplyr::ungroup() %>%
    # spread out annotation values
    tidyr::pivot_wider(
      id_cols = c(pair, score, score_unweighted),
      names_from = prediction_method,
      values_from = value
    ) %>%
    # get CSs and symbols
    dplyr::right_join(
      txv_master %>% dplyr::select(chrom, start = start.variant, end = end.variant,
                                   pair, variant, cs, enst, symbol),
      by = "pair"
    ) %>%
    # multiply score by expression binary
    dplyr::mutate(
      score_expressed = score * g_expression
    ) %>%
    # select columns
    dplyr::select(
      chrom, start, end,
      variant,
      cs,
      enst,
      symbol,
      score,
      where(is.numeric)
    )

  # predictions to save
  predictions <- scores %>%
    # filter to max score per CS-gene pair
    dplyr::group_by(cs, symbol) %>%
    dplyr::filter(score == max(score)) %>% # TODO: this will not get all unweighted maximums! must fix
    # annotate max gene per CS
    dplyr::group_by(cs) %>%
    dplyr::mutate(max = score == max(score)) %>%
    # order
    dplyr::select(chrom:end, variant, cs, symbol, score, max) %>%
    dplyr::arrange(-score)

  # max prediction per CS to save
  max <- predictions %>%
    dplyr::filter(max) %>%
    dplyr::select(variant, cs, symbol)

  # write tables
  write_tibble(scores, filename = out$Annotations)
  write_tibble(predictions, filename = out$Predictions)
  write_tibble(max, filename = out$MaxPredictions)
  }

  # 6) PERFORMANCE ======================================================================================================
  if(do_performance){
  cat("6) Measuring tgp performance...\n")
  # Generate PR curves (model performance metric)
  performance <- scores %>%
    # get performance
    get_PR(txv_master, drivers) %>%
    # add annotation level info
    purrr::map(~ dplyr::mutate(., level = sub("_.*", "", prediction_method)))

  pdf(out$PR, height = 10, width = 10, onefile = T)
  # PR score + max
  print(
    performance %>%
      plot_PR +
      ggplot2::ggtitle(out$Base %>% gsub("/", " ", .))
  )
  # PR max
  print(
    performance %>%
      purrr::map(~ .x %>% dplyr::filter(prediction_type == "max")) %>%
      plot_PR +
      ggplot2::ggtitle(out$Base %>% gsub("/", " ", .) %>% paste("max"))
  )
  # AUPRC
  print(
    performance %>%
      plot_AUPRC +
      ggplot2::ggtitle(out$Base %>% gsub("/", " ", .)) +
      ggplot2::coord_flip()
  )

  dev.off()

  # write table
  write_tibble(performance$summary, filename = out$Performance)
  }

  # 7) XGBoost MODEL TRAINING ======================================================================================================
  if(do_XGBoost){
  cat("7) Training an XGBoost model...\n")

  # format training set
  full <- scores %>%
    dplyr::mutate(label = (symbol %in% drivers$symbol) %>% as.numeric) %>%
    dplyr::group_by(cs) %>%
    dplyr::filter(any(label == 1)) %>%
    dplyr::ungroup()
  train <- list(data = full %>% dplyr::select(names(master)) %>% as.matrix,
                label = full$label)
  dtrain <- xgboost::xgb.DMatrix(data = train$data,
                                 label = train$label)
  # model training
  xgb1 <- xgboost::xgboost(data = dtrain,
                           max.depth = 5,
                           eta = 1,
                           nrounds = 500,
                           objective = "binary:logistic",
                           verbose = 1)
  # view feature importance plot
  xgb1_feature_importance_mat <- xgboost::xgb.importance(feature_names = colnames(dtrain), model = xgb1)
  xgb1_feature_importance_plot <- xgboost::xgb.ggplot.importance(importance_matrix = xgb1_feature_importance_mat)
  # save plot
  pdf(paste0(out$Base, "xgb_model_feature_importance.pdf"))
  print(xgb1_feature_importance_plot)
  dev.off()
  }

  # 8) SAVE ===
  # save(master,
  #      predictions,
  #      performance,
  #      xgb1,
  #      file = paste0(out$Base, "data.Rdata"))

  return(MA)
}



# 5) WEIGHTING ======================================================================================================
# # MANUAL WEIGHTING ===
# if(do_manual_weighting){
# celltype_of_interest <- unique(enriched$celltypes$name)
# if(length(celltype_of_interest) == 1){
#   to_add <- c("v_DHSs_signal",
#               "v_DHSs_specificity",
#               "txv_contact_ChIAPET_binary",
#               "txc_n_multicontact_binary_ChIAPET",
#               "gxv_specific_DHSs_closest_specific_genes",
#               "txv_exon")
#   to_multiply <- c("txv_TADs",
#                    "g_expression")
#   # celltype_of_interest = "BRST.HMEC" ; n_unique_manual_weights = 1
#   manual_models <- weight_and_score_manually(MA,
#                                              celltype_of_interest,
#                                              txv_master,
#                                              drivers,
#                                              to_add,
#                                              to_multiply,
#                                              n_unique_manual_weights)
#   write_tibble(manual_models, paste0(out$Base, "manual_weighting_models_performance.tsv"))
# } else { cat("dplyr::n_distinct(enriched$celltypes$name) != 1\nFunction will not work\n") }
# }


