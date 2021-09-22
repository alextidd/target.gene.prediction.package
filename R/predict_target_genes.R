#' Predict target genes of a list of fine-mapped non-coding variants for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param varfile A header-less BED file of fine-mapped trait-relevant variants, with metadata columns for the variant's name and the credible set it belongs to (chrom, start, stop, variant, cs)
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
predict_target_genes <- function(varfile,
                                 trait = NULL,
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
  dir.create(outDir)
  out <- list()
  out$Base <- paste0(outDir,"/") %>% { if(!is.null(trait)) paste0(., trait, "_") else . }
  out$Annotations <- paste0(out$Base, "target_gene_annotations.tsv")
  out$Predictions <- paste0(out$Base, "target_gene_predictions.tsv")
  out$Performance <- paste0(out$Base, "performance.tsv")
  out$PR <- paste0(out$Base, "PrecisionRecall.pdf")

  # import the variants
  cat("Importing variants...")
  variants <- target.gene.prediction.package::import_BED(
    variantsFile,
    metadata_cols = c("variant", "cs"))

  # import drivers and check that all symbols are in the GENCODE data
  cat("Importing driver genes...")
  drivers <- target.gene.prediction.package::read_tibble(driversFile)$V1
  target.gene.prediction.package::check_driver_symbols(drivers, driversFile)

  # import the contact data
  contact <- readRDS(paste0(referenceDir, "contact/contact.rda"))
  contact_metadata <- readRDS(paste0(referenceDir, "contact/contact_metadata.rda"))

  # import the DHSs data
  DHSs <- readRDS(paste0(referenceDir, "DHSs/DHSs.rda")) %>%
    target.gene.prediction.package::recursively_bind_rows(nest_names = c("Method", "Mark", "CellType"))
  DHSs_metadata <- readRDS(paste0(referenceDir, "DHSs/DHSs_metadata.rda"))

  # configure base weightings (not factoring in celltype enrichment)
  weight <- list(
    v = list(
      n_genes = 1),
    gxv = list(
      distance_score = 1,
      distance_rank = 1,
      closest = 1,
      contact = 2),
    gxc = list(
      multicontact = 3),
    gxd = list(
      multicontact = 3)
  )

  # ======================================================================================================
  # #### 1) CELL TYPE ENRICHMENT ####
  cat("1) Cell type enrichment...\n")
  enriched <- target.gene.prediction.package::get_enriched(DHSs, DHSs_metadata, variants)
  cat("Enriched cell types: ", enriched$celltypes$code, "\n")
  cat("Enriched tissues:", unique(enriched$celltypes$tissue), "\n")

  # Subset DHS annotations
  enriched_DHSs <- DHSs %>%
    dplyr::filter(CellType %in% enriched$tissues$code)
  enriched_contact_elements <- contact_metadata %>%
    dplyr::filter(Tissue %in% enriched$tissues$tissue) %>%
    dplyr::pull(list_element)

  # ======================================================================================================
  # #### 2) ENHANCER VARIANTS ####
  # get variants at DHSs ('enhancer variants')
  cat("1) Finding enhancer variants...\n")
  open_variants <- DHSs %>%
    dplyr::select(chrom:end) %>%
    dplyr::distinct() %>%
    target.gene.prediction.package::bed_intersect_left(
      variants, .,
      suffix = c("", ".DHS"),
      keepBcoords = T,
      keepBmetadata = F) %>%
    # unique DHS IDs to annotate
    dplyr::mutate(DHS = paste0(chrom, ":", start.DHS, "-", end.DHS))

  # ======================================================================================================
  # #### 3a) ENHANCER-LEVEL INPUTS ####
  # (Intersection with list of annotations, enhancers, GWAS statistics(?), etc...)
  v <- get_v_level_annotations()
  # bind and widen all variant-level annotations
  weighted_v_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = "variant",
    annotation.level = "v",
    v
  )

  # ======================================================================================================
  # #### 3b) GENE-LEVEL INPUTS ####
  g <- get_g_level_annotations()
  # bind and widen all gene-level annotations
  weighted_g_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = "enst",
    annotation.level = "g",
    g
  )

  # ======================================================================================================
  # #### 3c) DHS-LEVEL INPUTS ####
  d <- get_d_level_annotations()
  # bind and widen all DHS-level annotations
  weighted_d_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = "DHS",
    annotation.level = "d",
    d
  )

  # ======================================================================================================
  # #### 3d) CS-LEVEL INPUTS ####
  c <- get_c_level_annotations()
  # bind and widen all CS-level annotations
  weighted_c_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = "cs",
    annotation.level = "c",
    c
  )

  # ======================================================================================================
  # #### 3e) GENE-X-VARIANT-LEVEL INPUTS ####
  # (HiC interaction, distance, etc...)
  # The package will consider all genes as potential targets of a variant in CS's whose ranges are within the max distance
  gxv <- get_gxv_level_annotations()
  # bind and widen all gene-x-variant-level annotations
  weighted_gxv_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = c("variant", "cs", "enst"),
    annotation.level = "gxv",
    gxv
  )

  # ======================================================================================================
  # #### 3f) GENE-X-CS-LEVEL INPUTS ####
  gxc <- get_gxc_level_annotations()
  weighted_gxc_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = c("cs", "enst"),
    annotation.level = "gxc",
    gxc
  )

  # ======================================================================================================
  # #### 3g) GENE-X-DHS-LEVEL INPUTS ####
  gxd <- get_gxd_level_annotations()
  # bind and widen annotations
  weighted_gxd_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = c("DHS", "enst"),
    annotation.level = "gxd",
    gxd
  )

  # ======================================================================================================
  # #### 4) ALL INPUTS ####
  # Master variant-gene pair table
  # -> wide-format (one row per gene-variant pair, one column per annotation)
  # -> only variant-gene combinations with at least one pair annotation (e.g. HiChIP interaction, nearest or within 2Mb) are included
  # -> pair ID columns: | variant | enst |
  # -> annotation columns: | g_* | v_* | c_* | d_* | gxv_* | gxc_* | gxd_* |
  cat("3) Generating master table of gene x", trait, "variant pairs, with all annotation levels (genes, variants, gene-variant pairs, gene-credible set pairs)...\n")
  master <- open_variants %>%
    dplyr::select(variant, cs, DHS) %>%
    # Get annotations for all genes with TSSs within variant_to_gene_max_distance of every variant
    dplyr::right_join(  weighted_gxv_annotations      ) %>%
    dplyr::left_join(   weighted_gxc_annotations      ) %>%
    dplyr::left_join(   weighted_gxd_annotations      ) %>%
    dplyr::left_join(   weighted_v_annotations        ) %>%
    dplyr::left_join(   weighted_g_annotations        ) %>%
    dplyr::left_join(   weighted_d_annotations        ) %>%
    dplyr::left_join(   weighted_c_annotations        ) %>%
    # Replace all NAs in numeric annotation value columns with 0
    dplyr::mutate_if(is.numeric, tidyr::replace_na, replace = 0) %>%
    # Add gene symbols and logical driver column (T/F)
    dplyr::left_join(target.gene.prediction.package::TSSs %>% dplyr::select(enst, symbol)) %>%
    dplyr::mutate(driver = symbol %in% drivers)

  # get annotation columns (all numeric columns in the master table)
  annotation_cols <- master %>% dplyr::select(where(is.numeric)) %>% names()

  # Generate overall scores (rowwise sum function for now ; TODO: change!)
  predictions <- master %>%
    dplyr::mutate(score = rowSums(dplyr::across(annotation_cols))) %>%
    dplyr::select(cs, symbol, score, driver) %>%
    dplyr::distinct()

  # save output tables
  cat("Saving prediction files...\n")
  master %>% write_tibble(out$Annotations)
  predictions %>% write_tibble(out$Predictions)

  # ======================================================================================================
  # #### 5) PERFORMANCE ####

  # # Performance summaries: confusion matrix, p value, precision, recall, sensitivity, specificity, F score
  # # -> prediction                                      (positive = score > mean(score))
  # # -> max                                             (positive = group_by(cs) %>% score == max(score))
  # # -> all individual annotations converted to binary  (positive = score > 0)
  # performance <- dplyr::tibble()
  # for(col in c("max", "prediction")){ cat(col, "\n")
  #   curr.performance <- testable_predictions %>%
  #     target.gene.prediction.package::get_performance(., get(col)) %>%
  #     dplyr::mutate(Method = col)
  #   performance <- dplyr::bind_rows(performance, curr.performance)
  # }
  # for(col in annotation_cols){ cat(col, "\n")
  #   curr.performance <- testable_predictions_full %>%
  #     dplyr::mutate(score = dplyr::case_when(get(col) == 0 ~ FALSE, TRUE ~ TRUE)) %>%
  #     target.gene.prediction.package::get_performance(., score) %>%
  #     dplyr::mutate(Method = col)
  #   performance <- dplyr::bind_rows(performance, curr.performance)
  # }
  # performance <- performance %>% dplyr::select(Method, dplyr::everything())
  # # write performance table
  # performance %>% write_tibble(out$Performance)

  # # score distribution histograms
  # pdf("out/scores_distribution_per_annotation.pdf", width = 20, onefile = T)
  # for(annotation in annotation_cols){
  #   cat(annotation, "\n")
  #   print(target.gene.prediction.package::plot_scores_distribution(testable_predictions_full, annotation))
  # }
  # dev.off()

  # Generate PR curves (model performance metric)
  performance <- target.gene.prediction.package::get_PR(predictions, score)
  # PR of each individual annotation (columns of master)
  performance_all <- target.gene.prediction.package::get_PR(master, annotation_cols)
  # add annotation level info
  performance_all$PR <- performance_all$PR %>% dplyr::mutate(level = sub("_.*", "", prediction_type))

  pdf(outPR, height = 10, onefile = T)
    performance$PR %>% target.gene.prediction.package::plot_PR()
    performance$PR %>%
      dplyr::select(prediction_type, AUC) %>%
      dplyr::distinct() %>%
      ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_type, AUC),
                                   y = AUC)) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Predictor", y = "PR AUC")
    performance_all$PR %>%
      plot_PR(colour = prediction_method)
    performance_all$PR %>%
      dplyr::select(prediction_method, AUC) %>%
      dplyr::mutate(level = prediction_method %>% gsub("_.*", "", .)) %>%
      dplyr::distinct() %>%
      ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_type, AUC),
                          y = AUC,
                          fill = level)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_flip() +
      ggplot2::labs(x = "Predictor", y = "PR AUC")
  dev.off()

  # ======================================================================================================
  # #### 6) XGBoost model training? ####
  # XGBoost input
  master %>%
    dplyr::select(label = driver,
                  annotation_cols)


}




