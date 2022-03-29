#' Predict target genes of fine-mapped variants for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param trait Optional. The name of the trait of interest.
#' @param out_dir The output directory in which to save the predictions. Default is "./out/{trait}/{celltypes}/".
#' @param variants_file A BED file of trait-associated variants grouped by association signal, for example SNPs correlated with an index variant, or credible sets of fine-mapped variants
#' @param known_genes_file Optional. The file containing a list of trait known gene symbols. If do_performance is TRUE, must provide a known_genes_file.
#' @param reference_panels_dir The directory containing the external, accompanying reference panels data.
#' @param celltype_of_interest Optional. The celltype(s) of interest for the trait. Only annotations in these celltypes will be used to make predictions. Argument(s) must match the names of celltypes in the metadata. Make sure the celltype of interest has coverage across all annotations (TADs, HiChIP, expression, H3K27ac) in the metadata table.
#' @param tissue_of_interest Optional. The tissue(s) of interest for the trait. Only annotations in these tissues will be used to make predictions.  Argument(s) must match the names of tissues in the metadata.
#' @param variant_to_gene_max_distance The maximum absolute distance (bp) across which variant-gene pairs are considered. Default is 2Mb. The HiChIP data is also already filtered to 2Mb.
#' @param max_n_known_genes_per_CS In performance analysis, the maximum number of known genes within variant_to_gene_max_distance of the credible set.
#' @param min_proportion_of_variants_in_top_H3K27ac A threshold proportion of variants that reside in the specific H3K27ac-x-DHSs of a celltype for that celltype to be considered enriched. Default is 5\% (0.05).
#' @param do_all_celltypes If TRUE, the package will combine annotations across all available cell types, not just those with enriched enhancer variants. Default is FALSE.
#' @param do_all_celltypes_in_enriched_tissue If TRUE, the package will combine all annotations for the tissue of the enriched celltype(s), not just the specifically enriched celltype(s). Default is TRUE. Make sure the enriched celltype has coverage across all annotations (TADs, HiChIP, expression, H3K27ac) in the metadata table.
#' @param do_scoring If TRUE, runs the scoring chunk of the script, which combines all of the constituent MAE annotations into one score per transcript-variant pair. Default is FALSE.
#' @param do_performance If TRUE, runs the performance chunk of the script, which measures the performance of the score and each of its constituent annotations in predicting known genes as the targets of nearby variants. Default is FALSE.
#' @param do_XGBoost If TRUE, runs the XGBoost chunk of the script, which generates a model to predict the targets of variants from all available annotations and rates the importance of each annotation. Default is FALSE.
#' @param do_timestamp If TRUE, will save output into a subdirectory timestamped with the data/time of the run.
#' @param HiChIP If you are repeatedly running predict_target_genes, you can load the HiChIP object from the reference_panels_dir into the global environment and pass it to the function to prevent redundant re-loading each call to predict_target_genes.
#' @param H3K27ac If you are repeatedly running predict_target_genes, you can load the H3K27ac object from the reference_panels_dir into the global environment and pass it to the function to prevent redundant re-loading with each call to predict_target_genes.
#' @return A MultiAssayExperiment object with one assay object per annotation, one row per variant-transcript pair and one column per cell type (or 'value' if it is a non-cell-type-specific annotation).
#' @export
predict_target_genes <- function(trait = NULL,
                                 out_dir = NULL,
                                 variants_file = NULL,
                                 known_genes_file = NULL,
                                 reference_panels_dir = NULL,
                                 celltype_of_interest = NULL,
                                 tissue_of_interest = NULL,
                                 variant_to_gene_max_distance = 2e6,
                                 max_n_known_genes_per_CS = Inf,
                                 min_proportion_of_variants_in_top_H3K27ac = 0.05,
                                 do_all_celltypes = F,
                                 do_all_celltypes_in_enriched_tissue = T,
                                 do_scoring = T,
                                 do_performance = T,
                                 do_XGBoost = F,
                                 do_timestamp = F,
                                 HiChIP = NULL,
                                 H3K27ac = NULL) {
  # capture function arguments (do not run when testing internally)
  args <- as.list(environment())[names(as.list(environment())) %ni% c("HiChIP", "H3K27ac")]

  # for testing internally:
  # setwd("/working/lab_jonathb/alexandT/tgp") ; HiChIP = NULL ; H3K27ac = NULL ; celltype_of_interest = NULL ; tissue_of_interest = NULL ; trait="BC" ; out_dir = "out/" ; variants_file="/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/output/BC/variants/Michailidou2017/FM_variants.bed" ; known_genes_file = "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/traits/output/BC/known_genes.txt" ; reference_panels_dir = "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/output/" ; variant_to_gene_max_distance = 2e6 ; max_n_known_genes_per_CS = Inf ; min_proportion_of_variants_in_top_H3K27ac = 0.05 ; do_all_celltypes = F ; do_all_celltypes_in_enriched_tissue = T ; do_scoring = T ; do_performance = T ; do_XGBoost = T ; do_timestamp = F  ; library(devtools) ; load_all()
  # internally restore run environment:
  # # args <- dget("out/BC_Michailidou2017_FM_variants/enriched_celltypes/arguments_for_predict_target_genes.R") ; list2env(args, envir=.GlobalEnv)

  # SETUP ======================================================================================================

  # metadata for all annotations
  metadata <- read_tibble(paste0(reference_panels_dir, "metadata.tsv"), header = T)

  # check options
  {
    if(is.null(reference_panels_dir)){stop("Must provide the path to the accompanying reference panels directory!")}
    if(is.null(variants_file)){"Must provide the path to a file of trait variants!"}
    if (do_XGBoost) { do_scoring <- T }
    if (do_performance & is.null(known_genes_file)) { stop("do_performance = TRUE but no known_genes_file provided! Performance analysis requires a known_genes_file.") }
    if (!is.null(tissue_of_interest)) { if(tissue_of_interest %ni% metadata$tissue) {
      stop("Provided tissue_of_interest '", tissue_of_interest, "' is not represented in the available data. Must be one of...\n",
           paste(unique(metadata$tissue), collapse = ", ")) } }
    if (!is.null(celltype_of_interest)) {
      coi_annotations <- metadata %>% dplyr::filter(celltype == celltype_of_interest) %>% dplyr::pull(object) %>% unique
      if (celltype_of_interest %ni% metadata$celltype) {
        stop("Provided celltype_of_interest '", celltype_of_interest, "' is not represented in the available data. Must be one of...\n",
             paste(unique(metadata$celltype), collapse = ", ")) }
      if (length(setdiff(metadata$object, coi_annotations)) > 0){
        stop("Provided celltype_of_interest '", celltype_of_interest, "' does not have a full panel of annotations available.\nMissing annotation(s) = ",
             paste(setdiff(metadata$object, coi_annotations), collapse = ", "),
             "\nEither re-run at tissue-level (tissue_of_interest = ", unique(metadata[metadata$celltype == celltype_of_interest,]$tissue),
             ") or generate celltype-level annotations for '", celltype_of_interest, "' in the reference data.") } }
  }

  # define the output
  out <- list(
    base = "",
    tissue_enrichment = "tissue_fisher_enrichments.tsv",
    annotations = "target_gene_annotations.tsv",
    predictions = "target_gene_predictions_full.tsv",
    max_predictions = "target_gene_predictions_max.tsv",
    MA = "MA.rds",
    performance = "performance.tsv",
    PR = "PrecisionRecall.pdf",
    args = "arguments_for_predict_target_genes.R"
  )
  if(is.null(out_dir)){
    out <- paste0("out/", trait, "/") %>%
      {
        if(do_all_celltypes) paste0(., "all_celltypes/")
        else if(!is.null(tissue_of_interest)) paste0(., tissue_of_interest, "_tissue/")
        else if(!is.null(celltype_of_interest)) paste0(., celltype_of_interest, "_celltype/")
        else if(do_all_celltypes_in_enriched_tissue) paste0(., "enriched_tissues/")
        else paste0(., "enriched_celltypes/")
      } %>%
      { if(do_timestamp) paste0(., format(Sys.time(), "%Y%m%d_%H%M%S"), "/") else . } %>%
      { purrr::map(out, function(x) paste0(., x)) }
  } else { out <- out %>% purrr::map(function(x) paste0(out_dir, "/", x)) }
  dir.create(out$base, showWarnings = F, recursive = T)

  # write run arguments to output
  dput(args, file = out$args)

  # import the user-provided variants
  cat(" > Importing variants...\n")
  variants <- import_BED(
    variants_file,
    metadata_cols = c("variant", "cs"))

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
  DHSs <- readRDS(paste0(reference_panels_dir, "DHSs/DHSs.rds"))

  # import the TADs data
  cat(" > Importing TAD data...\n")
  TADs <- readRDS(paste0(reference_panels_dir, "TADs/TADs.rds"))

  # 1) CELL TYPE ENRICHMENT ======================================================================================================
  cat("1) Performing cell type enrichment...\n")
  enriched <- get_enriched(variants,
                           H3K27ac_specificity_ranked,
                           H3K27ac,
                           expression,
                           expressed,
                           HiChIP,
                           TADs,
                           metadata,
                           out,
                           min_proportion_of_variants_in_top_H3K27ac,
                           # options to manually change annotation groupings (passed by user):
                           celltype_of_interest,
                           tissue_of_interest,
                           do_all_celltypes,
                           do_all_celltypes_in_enriched_tissue)

  # 2) VARIANTS ======================================================================================================
  # get variant-to-gene universe ('enhancer variants')
  cat("2) Finding all genes near variants...\n")

  # The transcript-x-variant universe (masterlist of all possible transcript x variant pairs < variant_to_gene_max_distance apart)
  vxt_master <- variants %>%
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
                    pair = paste0(cs, ":", variant, "|", enst.TSS),
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
                               H3K27ac,
                               enriched,
                               vxt_master,
                               DHSs)

  # 3b) TRANSCRIPT-LEVEL INPUTS ====
  cat(" > T\tAnnotating transcripts...\n")
  t <- get_t_level_annotations(TSSs,
                               enriched)

  # 3c) GENE-LEVEL INPUTS ===
  cat(" > G\tAnnotating genes...\n")
  g <- get_g_level_annotations(vxt_master,
                               enriched)

  # 3d) CS-LEVEL INPUTS ====
  cat(" > C\tAnnotating credible sets...\n")
  c <- get_c_level_annotations(variants)

  # 3e) VARIANT-x-TRANSCRIPT-LEVEL INPUTS ====
  cat(" > VxT\tAnnotating variant x transcript pairs...\n")
  vxt <- get_vxt_level_annotations(variants,
                                   vxt_master,
                                   variant_to_gene_max_distance,
                                   enriched)

  # 3f) VARIANT-X-GENE-LEVEL INPUTS ===
  cat(" > VxG\tAnnotating variant x gene pairs...\n")
  vxg <- get_vxg_level_annotations(variants,
                                   vxt_master,
                                   enriched)

  # 3g) CS-X-TRANSCRIPT-LEVEL INPUTS ====
  cat(" > CxT\tAnnotating credible set x transcript pairs...\n")
  cxt <- get_cxt_level_annotations(vxt,
                                   variants)

  # 4) ALL INPUTS ======================================================================================================
  # Master variant-x-transcript matrix list
  # -> wide-format (one row per cs:variant|transcript pair, one column per celltype)
  # -> only variant-transcript combinations within 2Mb are included
  # -> pair ID rownames: variant|enst
  # Each annotation will be aggregated across samples, taking the maximum value per pair
  cat("4) Generating master table of transcript x", trait, "variant pairs, with all annotation levels...\n")
  master <- c(v, t, g, c, vxt, vxg, cxt) %>%
    purrr::map(~ matricise_by_pair(., vxt_master))

  # MultiAssayExperiment colData
  colData <- master %>%
    lapply(colnames) %>%
    unlist %>%
    as.data.frame %>%
    dplyr::distinct() %>%
    dplyr::rename(celltype = ".") %>%
    dplyr::left_join(metadata %>% dplyr::distinct(celltype, tissue), by = "celltype") %>%
    tibble::column_to_rownames("celltype")

  # celltype-specific annotations have as many columns as there are annotated celltypes
  # non-celltype-specific annotations have one `value` column, which applies across all celltypes
  MA <- MultiAssayExperiment::MultiAssayExperiment(experiments = master, colData = colData)
  saveRDS(MA, file = out$MA)

  # 5) SCORING ======================================================================================================
  if(do_scoring){
  cat("5) Scoring enhancer-gene pairs...\n")

  # declare weights (ALL OTHERS = 0.33)
  weights <- list(
    vxt_TADs = 1,
    vxt_inv_distance = 1,
    vxg_specific_H3K27ac_closest_specific_genes = 0, # 1
    vxt_intron = 1,
    vxg_missense = 1,
    vxg_nonsense = 1,
    vxg_splicesite = 1,
    t_H3K27ac_signal = 1,
    g_expressed = 1,
    g_expression_signal = 1,
    g_expression_specificity = 1,
    v_inv_n_genes = 0.66,
    c_inv_n_variants = 0.66
  )

  # Generating a single score for each variant-transcript pair, with evidence
  # mean(all non-celltype-specific values, enriched celltype-specific annotations) * (expression binary)
  scores <- master %>% names %>%
    sapply(function(annotation) {
      # weight annotations
      master[[annotation]] %>%
        tibble::as_tibble(rownames = "pair") %>%
        dplyr::mutate(
          value_unweighted = value,
          value = value * { if (annotation %in% names(weights)) weights[[annotation]] else 0.33 }
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
      vxt_master %>% dplyr::select(chrom, start = start.variant, end = end.variant,
                                   pair, variant, cs, enst, symbol),
      by = "pair"
    ) %>%
    # multiply score by expression binary
    dplyr::mutate(
      score_expressed = score * g_expressed
    ) %>%
    # select columns
    dplyr::select(
      cs,
      chrom:end,
      variant,
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
    dplyr::select(cs, chrom:end, variant, symbol, score, max) %>%
    dplyr::arrange(-score)

  # max prediction per CS to save
  max <- predictions %>%
    dplyr::filter(max) %>%
    dplyr::select(cs, variant, symbol)

  # write tables
  write_tibble(scores, filename = out$annotations)
  write_tibble(predictions, filename = out$predictions)
  write_tibble(max, filename = out$max_predictions)
  }

  # 6) PERFORMANCE ======================================================================================================
  if(do_performance){
  cat("6) Measuring tgp performance...\n")

  # import user-provided known genes and check that all symbols are in the GENCODE database
  cat(" > Importing known genes...\n")
  known_genes <- read_tibble(known_genes_file)$V1 %>%
    check_known_genes(known_genes_file)

  # Generate PR curves (model performance metric) (only testing protein-coding genes)
  performance <- scores %>%
    # get performance
    get_PR(vxt_master, known_genes, pcENSGs, max_n_known_genes_per_CS) %>%
    # add annotation level info
    purrr::map(~ dplyr::mutate(., level = sub("_.*", "", prediction_method)))

  pdf(out$PR, height = 10, width = 15, onefile = T)
    # PR score + max
    print(
      performance %>%
        plot_PR +
        ggplot2::ggtitle(out$base %>% gsub("/", " ", .))
    )
    # PR max
    print(
      performance %>%
        purrr::map(~ .x %>% dplyr::filter(prediction_type == "max")) %>%
        plot_PR +
        ggrepel::geom_text_repel(data = . %>%
                                   dplyr::ungroup() %>%
                                   dplyr::top_n(5, PR_AUC),
                                 min.segment.length = 0) +
        ggplot2::ggtitle(out$base %>% gsub("/", " ", .) %>% paste("max"))
    )
    # AUPRC
    print(
      performance %>%
        plot_AUPRC +
        ggplot2::ggtitle(out$base %>% gsub("/", " ", .))
    )
    print(
      performance$summary %>%
        dplyr::mutate(AUPRC = PR_AUC) %>%
        tidyr::pivot_longer(cols = c(Precision, Recall, AUPRC),
                            names_to = "metric",
                            values_to = "performance") %>%
        ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_method, PR_AUC),
                                     y = performance,
                                     fill = level)) +
        ggplot2::geom_col() +
        ggplot2::facet_wrap(~metric) +
        ggplot2::coord_flip() +
        ggplot2::theme_classic() +
        ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
        ggsci::scale_fill_igv()
    )
  dev.off()

  # write table
  write_tibble(performance$summary, filename = out$performance)
  }

  # 7) XGBoost MODEL TRAINING ======================================================================================================
  if(do_XGBoost){
  cat("7) Training an XGBoost model...\n")

  # format training set
  full <- scores %>%
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
                           nrounds = 500,
                           objective = "binary:logistic",
                           verbose = 1)
  # view feature importance plot
  xgb1_feature_importance_mat <- xgboost::xgb.importance(feature_names = colnames(dtrain), model = xgb1)
  xgb1_feature_importance_plot <- xgboost::xgb.ggplot.importance(importance_matrix = xgb1_feature_importance_mat)
  # save plot
  pdf(paste0(out$base, "xgb_model_feature_importance.pdf"))
  print(xgb1_feature_importance_plot)
  dev.off()
  }

  # 8) SAVE ===
  # save(master,
  #      predictions,
  #      performance,
  #      xgb1,
  #      file = paste0(out$base, "data.Rdata"))

  return(MA)
}

