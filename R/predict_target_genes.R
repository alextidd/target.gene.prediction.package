#' Predict target genes of a list of fine-mapped non-coding variants for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param trait The name of the trait of interest
#' @param tissue The tissue(s) of action for the trait
#' @param outDir The output directory in which to save the predictions. Default is current directory
#' @param referenceDir The directory containing the external, accompanying reference data
#' @param variantsFile A BED file of trait-associated variants grouped by association signal, for example SNPs correlated with an index variant, or credible sets of fine-mapped variants
#' @param driversFile The file containing a list of driver gene symbols
#' @param variant_to_gene_max_distance The maximum absolute distance (bp) across which variant-gene pairs are considered. Default is 2Mb.
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

  # silence "no visible binding" NOTE for data variables
  . <- variant <- enst <- start.variant <- start.TSS <- variant.variant <- enst.TSS <- score <- annotation <- estimate <-
    p.value <- enst.query <- annotation.annotation <- variant.query <- n <- end <- pair.score <-
    NULL

  # define the output
  dir.create(outDir, recursive = T, showWarnings = F)
  out <- list()
  out$Base <- paste0(outDir,"/") %>% { if(!is.null(trait)) paste0(., trait, "_") else . }
  out$Annotations <- paste0(out$Base, "target_gene_annotations.tsv")
  out$Predictions <- paste0(out$Base, "target_gene_predictions.tsv")
  out$Performance <- paste0(out$Base, "performance.tsv")
  out$PR <- paste0(out$Base, "PrecisionRecall.pdf")

  # import the variants
  cat("Importing variants...\n")
  variants <- target.gene.prediction.package::import_BED(
    variantsFile,
    metadata_cols = c("variant", "cs"))

  # import drivers and check that all symbols are in the GENCODE data
  cat("Importing driver genes...\n")
  drivers <- target.gene.prediction.package::read_tibble(driversFile)$V1
  target.gene.prediction.package::check_driver_symbols(drivers, driversFile)

  # import the contact data
  contact <- readRDS(paste0(referenceDir, "contact/contact.rda"))
  contact_metadata <- readRDS(paste0(referenceDir, "contact/contact_metadata.rda"))

  # import the DHSs data
  DHSs <- readRDS(paste0(referenceDir, "DHSs/DHSs.rda"))
  DHSs_metadata <- readRDS(paste0(referenceDir, "DHSs/DHSs_metadata.rda"))
  DHSs_master <-  DHSs[[1]] %>%
    dplyr::distinct(chrom, start, end, DHS) %>%
    dplyr::mutate(DHSID = paste0(".", dplyr::row_number()))

  # import the TADs data
  TADs <- readRDS(paste0(referenceDir, "TADs/TADs.rda"))

  # configure base weightings (not factoring in celltype enrichment)
  weight <- list(
    v = list(
      n_genes = 1),
    gxv = list(
      distance_score = 1,
      distance_rank = 1,
      closest = 1,
      contact = 2,
      promoter = 5,
      exon = 5,
      intron = 5),
    gxc = list(
      multicontact = 3),
    gxd = list(
      multicontact = 3,
      closest = 1)
  )

  # ======================================================================================================
  # #### 1) CELL TYPE ENRICHMENT ####
  cat("1) Cell type enrichment...\n")
  enriched <- target.gene.prediction.package::get_enriched(DHSs$specificity, DHSs_metadata, variants)
  cat("Enriched cell type(s): ", enriched$celltypes$code, "\n")
  cat("Enriched tissue(s):", unique(enriched$celltypes$tissue), "\n")

  # Subset DHS annotations
  enriched$DHSs <- DHSs %>%
    lapply(dplyr::select, chrom:DHS, dplyr::contains(enriched$tissues$code))
  enriched$contact_elements <- contact_metadata %>%
    dplyr::filter(Tissue %in% enriched$tissues$tissue) %>%
    dplyr::pull(list_element)

  # ======================================================================================================
  # #### 2) ENHANCER VARIANTS ####
  # get variants at DHSs ('enhancer variants')
  cat("2) Finding enhancer variants...\n")
  open_variants <- DHSs_master %>%
    target.gene.prediction.package::bed_intersect_left(
      variants, .,
      keepBcoords = F,
      keepBmetadata = T)

  # The SNP-gene universe (masterlist of all possible gene x open variant pairs <2Mb apart)
  gxv_master <- open_variants %>%
    valr::bed_slop(both = variant_to_gene_max_distance,
                   genome = target.gene.prediction.package::ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs,
                        suffix = c(".variant", ".TSS")) %>%
    dplyr::distinct(variant = variant.variant,
                    cs = cs.variant,
                    DHS = DHS.variant,
                    enst = enst.TSS
                    ) %>%
    dplyr::mutate(pair = paste0(variant, "|", enst))

  # ======================================================================================================
  cat("3) Scoring enhancer-gene pairs...\n")
  # #### 3a) ENHANCER-LEVEL INPUTS ####
  # (Intersection with list of annotations, enhancers, GWAS statistics(?), etc...)
  v <- get_v_level_annotations(open.variants = open_variants,
                               enriched.DHSs = enriched$DHSs,
                               variant.to.gene.max.distance = variant_to_gene_max_distance)

  # ======================================================================================================
  # #### 3b) GENE-LEVEL INPUTS ####
  g <- get_g_level_annotations(enriched.DHSs = enriched$DHSs)

  # ======================================================================================================
  # #### 3c) DHS-LEVEL INPUTS ####
  d <- get_d_level_annotations(open.variants = open_variants)

  # ======================================================================================================
  # #### 3d) CS-LEVEL INPUTS ####
  c <- get_c_level_annotations(open.variants = open_variants)

  # ======================================================================================================
  # #### 3e) GENE-X-VARIANT-LEVEL INPUTS ####
  # (3C interaction, distance, etc...)
  # The package will consider all genes as potential targets of a variant in CS's whose ranges are within the max distance
  gxv <- get_gxv_level_annotations(open.variants = open_variants,
                                   variant.to.gene.max.distance = variant_to_gene_max_distance,
                                   .enriched = enriched,
                                   .contact = contact)

  # ======================================================================================================
  # #### 3f) GENE-X-CS-LEVEL INPUTS ####
  gxc <- get_gxc_level_annotations(.gxv = gxv,
                                   open.variants = open_variants)

  # ======================================================================================================
  # #### 3g) GENE-X-DHS-LEVEL INPUTS ####
  gxd <- get_gxd_level_annotations(ensts.near.vars = unique(gxv$gxv_inv_distance$enst),
                                   DHSs.master = DHSs_master)

  # ======================================================================================================
  # #### 4) ALL INPUTS ####
  # Master variant-gene pair table
  # -> wide-format (one row per gene-variant pair, one column per annotation)
  # -> only variant-gene combinations with at least one pair annotation (e.g. HiChIP interaction, nearest or within 2Mb) are included
  # -> pair ID columns: | variant | enst |
  # -> annotation columns: | g_* | v_* | c_* | d_* | gxv_* | gxc_* | gxd_* |
  cat("3) Generating master table of gene x", trait, "variant pairs, with all annotation levels...\n")
  master <- c(v, g, d, c, gxv, gxc, gxd) %>%
    purrr::map(~ matricise_by_pair(., gxv_master))
  master %>% write_tibble(out$Annotations)

  # MultiAssayExperiment
  samplelist <- read.delim("/working/lab_georgiat/jonathB/PROJECTS/GitHub_repos/target.gene.prediction.package/JB/coldata.txt", stringsAsFactors = F)
  col_from_assays <- lapply(master, colnames) %>% unlist() %>% as.data.frame() %>% dplyr::distinct()
  coldata <- dplyr::left_join(col_from_assays, samplelist, by = c("."="X"))
  rownames(coldata) <- coldata$.

  MA <- MultiAssayExperiment::MultiAssayExperiment(experiments = master, colData = coldata)
  return(MA)

}



# open_variants %>%
#   dplyr::select(variant, cs, DHS) %>%
#   # Get annotations for all genes with TSSs within variant_to_gene_max_distance of every variant
#   dplyr::right_join( weighted_gxv_annotations, by = c("variant", "cs")) %>%
#   dplyr::left_join(  weighted_gxc_annotations, by = c("cs", "enst")) %>%
#   dplyr::left_join(  weighted_gxd_annotations, by = c("DHS", "enst")) %>%
#   dplyr::left_join(  weighted_v_annotations,   by = "variant") %>%
#   dplyr::left_join(  weighted_g_annotations,   by = "enst") %>%
#   dplyr::left_join(  weighted_d_annotations,   by = "DHS") %>%
#   dplyr::left_join(  weighted_c_annotations,   by = "cs") %>%
#   # Replace all NAs in numeric annotation value columns with 0
#   dplyr::mutate_if(is.numeric, tidyr::replace_na, replace = 0) %>%
#   # Add gene symbols and logical driver column (T/F)
#   dplyr::left_join(target.gene.prediction.package::TSSs %>% dplyr::select(enst, symbol)) %>%
#   dplyr::mutate(driver = symbol %in% drivers)
#
# # get annotation columns (all numeric columns in the master table)
# annotation_cols <- master %>% dplyr::select(where(is.numeric)) %>% names()
#
# # Generate overall scores (rowwise sum function for now ; ## TODO: change!)
# predictions <- master %>%
#   dplyr::mutate(score = rowSums(dplyr::across(dplyr::all_of(annotation_cols)))) %>%
#   dplyr::select(cs, symbol, score, driver) %>%
#   dplyr::distinct()
# predictions_out <- predictions %>%
#   dplyr::select(-driver) %>%
#   dplyr::group_by(cs, symbol) %>%
#   dplyr::filter(score == max(score))
#
# # save output tables
# cat("Saving prediction files...\n")
# master %>% write_tibble(out$Annotations)
# predictions_out %>% write_tibble(out$Predictions)
#
# # ======================================================================================================
# # #### 5) PERFORMANCE ####
#
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
#
# # ======================================================================================================
# # #### 6) XGBoost model training? ####
# # XGBoost input
# train <- master %>%
#   dplyr::rename(label = driver)
#
# # MAE reformatting
# # 1. Intersect DHS masterlist with variants
# mat_var_DHSs <- open_variants %>%
#   dplyr::distinct(variant, DHSID, DHS) %>%
#   dplyr::left_join(., gxv_master)

