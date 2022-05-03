#' Predict target genes of fine-mapped variants for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param trait Optional. The name of the trait of interest.
#' @param out_dir Optional. The output directory in which to save the predictions. Default is "./out/{trait}/{celltypes}/".
#' @param sub_dir Optional. The output subdirectory to be appended to the default path ("./out/{trait}/{celltypes}/{sub_dir}/").
#' @param variants_file A BED file of trait-associated variants grouped by association signal, for example SNPs correlated with an index variant, or credible sets of fine-mapped variants
#' @param known_genes_file Optional. The file containing a list of trait known gene symbols. If do_performance is TRUE, must provide a known_genes_file.
#' @param reference_panels_dir The directory containing the external, accompanying reference panels data.
#' @param celltype_of_interest Optional. The celltype(s) of interest for the trait. Only annotations in these celltypes will be used to make predictions. Argument(s) must match the names of celltypes in the metadata. Make sure the celltype of interest has coverage across all annotations (TADs, HiChIP, expression, H3K27ac) in the metadata table.
#' @param tissue_of_interest Optional. The tissue(s) of interest for the trait. Only annotations in these tissues will be used to make predictions.  Argument(s) must match the names of tissues in the metadata.
#' @param variant_to_gene_max_distance The maximum absolute distance (bp) across which variant-gene pairs are considered. Measured as the distance between the variant and the gene's TSS. Default is 2Mb. The HiChIP data is also already filtered to 2Mb.
#' @param max_n_known_genes_per_CS In performance analysis, the maximum number of known genes within variant_to_gene_max_distance of the credible set.
#' @param celltypes Dictates which celltypes' annotations are used. Must be one of c("enriched_celltypes", "enriched_tissues", "all_celltypes"). If "enriched_celltypes", annotations from only the enriched celltype(s) will be used. The enriched celltype(s) must have coverage across all annotations (TADs, HiChIP, expression, H3K27ac) in the metadata table for this to work. If "enriched_tissues", all annotations from the tissue of the enriched celltype(s) will be used. If "all_celltypes", the enrichment analysis is skipped and annotations from all available celltypes will be used. Default is "enriched_tissues".
#' @param do_scoring If TRUE, runs the scoring chunk of the script, which combines all of the constituent MAE annotations into one score per transcript-variant pair. Default is FALSE.
#' @param do_performance If TRUE, runs the performance chunk of the script, which measures the performance of the score and each of its constituent annotations in predicting known genes as the targets of nearby variants. Default is FALSE.
#' @param do_XGBoost If TRUE, runs the XGBoost chunk of the script, which generates a model to predict the targets of variants from all available annotations and rates the importance of each annotation. Default is FALSE.
#' @param do_timestamp If TRUE, will save output into a subdirectory timestamped with the data/time of the run.
#' @param HiChIP If you are repeatedly running predict_target_genes, you can load the HiChIP object from the reference_panels_dir into the global environment and pass it to the function to prevent redundant re-loading each call to predict_target_genes.
#' @param H3K27ac If you are repeatedly running predict_target_genes, you can load the H3K27ac object from the reference_panels_dir into the global environment and pass it to the function to prevent redundant re-loading with each call to predict_target_genes.
#' @return The `annotations` tibble with one row per variant-x-transcript pair, one column per annotation and the resulting weighted score for each pair.
#' @export
predict_target_genes <- function(trait = NULL,
                                 out_dir = NULL,
                                 sub_dir = NULL,
                                 variants_file = NULL,
                                 known_genes_file = NULL,
                                 reference_panels_dir = NULL,
                                 celltype_of_interest = NULL,
                                 tissue_of_interest = NULL,
                                 celltypes = "enriched_tissues",
                                 variant_to_gene_max_distance = 2e6,
                                 max_n_known_genes_per_CS = Inf,
                                 do_scoring = T,
                                 do_performance = T,
                                 do_XGBoost = F,
                                 do_timestamp = F,
                                 HiChIP = NULL,
                                 H3K27ac = NULL) {

  # capture function arguments (do not run when testing internally)
  args <- as.list(environment())[names(as.list(environment())) %ni% c("HiChIP", "H3K27ac")]

  # for testing internally:
  # setwd("/working/lab_jonathb/alexandT/tgp") ; trait="BC_Michailidou2017_FM" ; celltypes = "enriched_tissues" ; variants_file=paste0("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/output/",trait,"/variants.bed") ; known_genes_file = paste0("/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/output/",trait,"/known_genes.txt") ; reference_panels_dir = "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/output/" ; variant_to_gene_max_distance = 2e6 ; max_n_known_genes_per_CS = Inf ; HiChIP = NULL ; H3K27ac = NULL ; celltype_of_interest = NULL ; tissue_of_interest = NULL ; out_dir = NULL ; sub_dir = NULL ; do_scoring = T ; do_performance = T ; do_XGBoost = T ; do_timestamp = F  ; library(devtools) ; load_all()
  # for internally restoring a previous run environment:
  # args <- dget("out/BC_Michailidou2017_FM/enriched_tissues/20220503_100027/arguments_for_predict_target_genes.R") ; list2env(args, envir=.GlobalEnv) ; library(devtools) ; load_all()

  # SETUP ======================================================================================================

  # metadata for all annotations
  metadata <- read_tibble(paste0(reference_panels_dir, "metadata.tsv"), header = T)

  # check arguments
  check_arguments(metadata,
                  variants_file,
                  known_genes_file,
                  reference_panels_dir,
                  celltype_of_interest,
                  tissue_of_interest,
                  celltypes,
                  do_scoring,
                  do_performance,
                  do_XGBoost)

  # define the output files
  enriched_dir <- {
    if (!is.null(celltype_of_interest)) paste0(celltype_of_interest, "_celltype")
    else if (!is.null(tissue_of_interest)) paste0(tissue_of_interest, "_tissue")
    else paste0(celltypes)
  }
  out <- list(
    base = "",
    tissue_enrichments = "tissue_enrichments.tsv",
    master = "master.rds",
    annotations = "target_gene_annotations.tsv",
    predictions_full = "target_gene_predictions_full.tsv",
    predictions_max = "target_gene_predictions_max.tsv",
    performance = "performance.tsv",
    PR_plot = "precision_recall_plot.pdf",
    args = "arguments_for_predict_target_genes.R",
    XGBoost = "XGBoost_feature_importance.tsv",
    XGBoost_plot = "XGBoost_feature_importance_plot.pdf"
  ) ; if (is.null(out_dir)) {
    out <- "out/" %>%
      paste0(trait, "/", enriched_dir, "/") %>%
      { if(do_timestamp) paste0(., format(Sys.time(), "%Y%m%d_%H%M%S"), "/") else . } %>%
      { if(!is.null(sub_dir)) paste0(., sub_dir, "/") else . } %>%
      { purrr::map(out, function(x) paste0(., x)) }
  } else { out <- out %>% purrr::map(function(x) paste0(out_dir, "/", x)) }
  dir.create(out$base, showWarnings = F, recursive = T)

  # write run arguments to output
  dput(args, file = out$args)

  # ggplot theme settings
  ggplot2::theme_set(
    ggplot2::theme_classic() +
      ggplot2::theme(
        title = ggplot2::element_text(size = 20),
        text  = ggplot2::element_text(size = 18)))

  # import the user-provided variants
  cat(" > Importing variants...\n")
  variants <- import_BED(
    variants_file,
    metadata_cols = c("variant", "cs")) %>%
    dplyr::mutate(variant = as.character(variant),
                  cs = as.character(cs))

  # import the HiChIP data
  if (is.null(HiChIP)) {
    cat(" > Importing HiChIP data...\n")
    HiChIP <- readRDS(paste0(reference_panels_dir, "HiChIP/HiChIP.rds"))
  }

  # import the H3K27ac-x-DHS binning data
  if (is.null(H3K27ac)) {
    cat(" > Importing H3K27ac-x-DHS binning data...\n")
    H3K27ac <- readRDS(paste0(reference_panels_dir, "H3K27ac/H3K27ac.rds"))
  }
  H3K27ac_specificity_ranked <- readRDS(paste0(reference_panels_dir, "H3K27ac/H3K27ac_specificity_rank.rds"))

  # import the expression data
  cat(" > Importing RNA-seq expression data...\n")
  expression <- readRDS(paste0(reference_panels_dir, "expression/expression.rds"))
  expressed <- readRDS(paste0(reference_panels_dir, "expression/expressed.rds"))

  # import DHSs master
  cat(" > Importing DHS data...\n")
  DHSs <- readRDS(paste0(reference_panels_dir, "DHSs/DHSs.rds"))

  # import the TADs data
  cat(" > Importing TAD data...\n")
  TADs <- readRDS(paste0(reference_panels_dir, "TADs/TADs.rds"))

  # 1) CELL TYPE ENRICHMENT ======================================================================================================
  cat("1) Performing cell type enrichment...\n")
  enriched <- get_enriched(variants,
                           DHSs,
                           H3K27ac_specificity_ranked,
                           H3K27ac,
                           expression,
                           expressed,
                           HiChIP,
                           TADs,
                           metadata,
                           out,
                           # options to manually choose annotation group (passed by user):
                           celltype_of_interest,
                           tissue_of_interest,
                           celltypes)

  # 2) VARIANTS ======================================================================================================
  # get variant-to-gene universe ('enhancer variants')
  cat("2) Finding all genes near variants...\n")

  # The transcript-x-variant universe (masterlist of all possible transcript x variant pairs < variant_to_gene_max_distance apart)
  vxt_master <- get_vxt_master(variants,
                               TSSs,
                               variant_to_gene_max_distance)

  # 3) ANNOTATING ======================================================================================================
  cat("3) Annotating variant-transcript pairs at every level...\n")

  cat("> V\tAnnotating variants...\n")
  v <- get_v_level_annotations(variants,
                               H3K27ac,
                               enriched,
                               vxt_master,
                               DHSs)

  cat(" > T\tAnnotating transcripts...\n")
  t <- get_t_level_annotations(TSSs,
                               DHSs,
                               enriched)

  cat(" > G\tAnnotating genes...\n")
  g <- get_g_level_annotations(vxt_master,
                               enriched)

  cat(" > C\tAnnotating credible sets...\n")
  c <- get_c_level_annotations(variants)

  cat(" > VxT\tAnnotating variant x transcript pairs...\n")
  vxt <- get_vxt_level_annotations(variants,
                                   DHSs,
                                   vxt_master,
                                   variant_to_gene_max_distance,
                                   enriched)

  cat(" > VxG\tAnnotating variant x gene pairs...\n")
  vxg <- get_vxg_level_annotations(variants,
                                   vxt_master,
                                   enriched)

  cat(" > CxT\tAnnotating credible set x transcript pairs...\n")
  cxt <- get_cxt_level_annotations(vxt,
                                   variants)

  cat(" > CxG\tAnnotating credible set x gene pairs...\n")
  cxg <- get_cxg_level_annotations(vxt_master)

  # 4) ALL INPUTS ======================================================================================================
  # Master variant-x-transcript matrix list
  # -> wide-format (one row per cs:variant|transcript pair, one column per celltype)
  # -> only variant-transcript combinations within 2Mb are included
  # -> pair ID rownames: variant|enst
  # Each annotation will be aggregated across samples, taking the maximum value per pair
  cat("4) Generating master table of transcript x", trait, "variant pairs, with all annotation levels...\n")
  master <- c(v, t, g, c, vxt, vxg, cxt, cxg) %>%
    purrr::map(~ matricise_by_pair(., vxt_master))

  # 5) SCORING ======================================================================================================
  if(do_scoring){
  cat("5) Scoring enhancer-gene pairs...\n")

  # get weights (ALL OTHERS = 0.33)
  annotations_metadata <- read_tibble("data/metadata.tsv", header = T)
  weights <- as.list(annotations_metadata$weight) %>% setNames(annotations_metadata$annotation)
  if(length(setdiff(names(master), names(weights))) > 0){stop("Annotation ", setdiff(names(master), names(weights)), " does not have a weight in the annotations metadata!")}

  # weight and average
  means <- names(master) %>%
    sapply(function(annotation){rowMeans(master[[annotation]] * weights[[annotation]])})
  means <- cbind(means,
                 score = rowMeans(means),
                 score_expressed = rowMeans(means) * means[, "g_expressed"])

  # predictions to save
  # all vxt annotations
  annotations <- dplyr::bind_cols(
    vxt_master %>% dplyr::select(cs, chrom, start = start.variant, end = end.variant, variant, symbol, enst),
    means %>% dplyr::as_tibble()) %>%
    dplyr::arrange(-score_expressed)
  # max score per vxg
  predictions_full <- annotations %>%
    dplyr::select(cs:symbol, score = score_expressed) %>%
    dplyr::group_by(variant, symbol) %>%
    dplyr::filter(score == max(score))
  # max gene per cs
  predictions_max <- predictions_full %>%
    dplyr::select(cs, variant, symbol, score) %>%
    dplyr::group_by(cs) %>%
    dplyr::filter(score == max(score) & score > 0)

  # write tables
  write_tibble(annotations, filename = out$annotations)
  write_tibble(predictions_full, filename = out$predictions_full)
  write_tibble(predictions_max, filename = out$predictions_max)
  }

  # 6) PERFORMANCE ======================================================================================================
  if(do_performance){
  cat("6) Measuring tgp performance...\n")

  # import user-provided known genes and check that all symbols are in the GENCODE database
  cat(" > Importing known genes...\n")
  known_genes <- read_tibble(known_genes_file)$V1 %>%
    check_known_genes(known_genes_file)

  # Generate PR curves (model performance metric) (only testing protein-coding genes)
  performance <- annotations %>%
    # get performance
    get_PR(vxt_master, known_genes, pcENSGs, max_n_known_genes_per_CS) %>%
    # add annotation level info
    purrr::map(~ dplyr::mutate(., level = sub("_.*", "", prediction_method)))

  # plot extras
  weight_facets <- dplyr::left_join(
    performance$summary %>% dplyr::select(prediction_method),
      dplyr::tibble(prediction_method = names(weights),
                    weight = unlist(weights)),
    by = "prediction_method") %>%
    dplyr::mutate(
      weight = factor(dplyr::case_when(
        grepl("^score", prediction_method) ~ "score",
        is.na(weight) ~ "0.33",
        TRUE ~ as.character(weight)),
        levels = c("score", "1", "0.66", "0.33")))

  title_plot <- list(ggplot2::labs(
    title = paste0(
      "\nTrait = ", trait,
      "\nmax n known genes per CS = ", max_n_known_genes_per_CS,
      "; max distance = ", variant_to_gene_max_distance,
      "\nEnrichment = ", enriched_dir
    ),
    subtitle = paste0(
      "Tissue(s) = ", enriched$celltypes$tissue %>% unique %>% paste(collapse = ", "),
      "\n",
      paste(strwrap(
        paste0(
          "Celltype(s) = ", enriched$celltypes$celltype %>% unique %>% paste(collapse = ", ")
        ),
        width = 70
      ), collapse = "\n")
    )
  ))

  {pdf(out$PR_plot, height = 10, width = 15, onefile = T)
    # PR score + max
    print(plot_PR(performance, colour = prediction_method) +
            title_plot)

    # PR max
    print(
      performance %>%
        purrr::map(~ .x %>% dplyr::filter(prediction_type == "max")) %>%
        plot_PR(colour = prediction_method) +
        ggrepel::geom_text_repel(min.segment.length = 0, max.overlaps = Inf) +
        title_plot
    )
    # F score max
    print(
      performance$summary %>%
        dplyr::filter(prediction_type == "max") %>%
        dplyr::mutate(level = prediction_method %>% gsub("_.*", "", .),
                      F_score = F_score %>% tidyr::replace_na(0)) %>%
        dplyr::distinct() %>%
        ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_method, F_score),
                                     y = F_score,
                                     fill = level)) +
        ggplot2::geom_col() +
        ggplot2::labs(x = "Predictor",
                      y = "F score") +
        ggsci::scale_fill_igv() +
        ggplot2::coord_flip() +
        title_plot
    )
    print(
      performance$summary %>%
        dplyr::left_join(weight_facets, by = "prediction_method") %>%
        dplyr::mutate(dplyr::across(where(is.numeric), tidyr::replace_na, 0),
                      fsc = F_score) %>%
        tidyr::pivot_longer(cols = c(F_score, PR_AUC, Precision, Recall),
                            names_to = "metric",
                            values_to = "performance") %>%
        ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_method, fsc),
                                     y = performance,
                                     fill = level)) +
        ggplot2::geom_col() +
        ggplot2::facet_grid(weight ~ metric,
                            scales = "free_y", space = "free_y") +
        ggplot2::coord_flip() +
        ggplot2::theme(axis.title = ggplot2::element_blank()) +
        ggsci::scale_fill_igv() +
        title_plot
    )
  dev.off()}

  # write table
  write_tibble(performance$summary, filename = out$performance)
  }

  # 7) XGBoost MODEL TRAINING ======================================================================================================
  if(do_XGBoost){
  cat("7) Training an XGBoost model...\n")

  # format training set
  full <- annotations %>%
    dplyr::mutate(label = (symbol %in% known_genes$symbol) %>% as.numeric) %>%
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
                           nthread = 1,
                           nrounds = 100,
                           objective = "binary:logistic",
                           verbose = 1)
  # view feature importance plot
  xgb1_feature_importance_mat <- xgboost::xgb.importance(feature_names = colnames(dtrain), model = xgb1)
  xgb1_feature_importance_plot <- xgboost::xgb.ggplot.importance(importance_matrix = xgb1_feature_importance_mat)
  # save plot
  {
    pdf(out$XGBoost_plot, width = 15, height = 15)
    print(xgb1_feature_importance_plot + title_plot)
    dev.off()
  }
  # write table
  write_tibble(xgb1_feature_importance_mat, filename = out$XGBoost)
  }

  cat("\nDONE!\n")
  # 8) SAVE ===
  return(annotations)
}

