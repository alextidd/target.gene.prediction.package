#' Predict target genes of a list of fine-mapped non-coding variants for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param varfile A header-less BED file of fine-mapped trait-relevant variants, with metadata columns for the variant's name and the credible set it belongs to (chrom, start, stop, variant, cs)
#' @param trait The name of the trait of interest.
#' @param tissue The tissue(s) of action for the trait.
#' @param outDir The output directory in which to save the predictions. Default is current directory.
#' @param referenceDir The directory containing the external, accompanying reference data.
#' @param driversFile The file containing a list of driver gene symbols
#' @param variant_to_gene_max_distance The maximum absolute distance (bp) across which variant-gene pairs are considered. Default is 2Mb.
#'
#' @return A file of variant-gene pair predictions, with associated scores, saved in the given output directory.
#' @export
predict_target_genes <- function(varfile,
                                 trait = NULL,
                                 tissue = NULL,
                                 outDir = "out",
                                 referenceDir = "/working/lab_georgiat/alexandT/target.gene.prediction.package/external_data/reference/",
                                 driversFile = "/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/Traits/BC/KnownGenes/breast_cancer_drivers_2021.txt",
                                 variant_to_gene_max_distance = 2e6){

  # for testing:
  # library(devtools) ; load_all() ; variants=BCVars %>% dplyr::select(chrom:cs); trait="BC"; outDir="test_output"; referenceDir="/working/lab_georgiat/alexandT/target.gene.prediction.package/external_data/reference/" ; driversFile = "/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/Traits/BC/KnownGenes/breast_cancer_drivers_2021.txt" ; variant_to_gene_max_distance=2e6

  # silence "no visible binding" NOTE for data variables
  . <- variant <- enst <- start.variant <- start.TSS <- variant.variant <- enst.TSS <- score <- annotation <- estimate <-
    p.value <- enst.query <- annotation.annotation <- variant.query <- n <- end <- pair.score <-
    NULL

  # define the output
  dir.create(outDir)
  outBase <- paste0(outDir,"/") %>% { if(!is.null(trait)) paste0(., trait, "_") else . }
  outPredictions <- paste0(outBase, "target_gene_predictions.tsv")
  outPR <- paste0(outBase, "PrecisionRecall.pdf")

  # import the variants
  variants <- target.gene.prediction.package::import_BED(
    varfile,
    metadata_cols = c("variant", "cs"))

  # import the contact data
  contact <- readRDS(paste0(referenceDir, "contact/contact.rda"))

  # import the DHSs data
  DHSs <- readRDS(paste0(referenceDir, "DHSs/DHSs.rda")) %>%
    target.gene.prediction.package::recursively_bind_rows(nest_names = c("Method", "Mark", "CellType"))

  # load drivers and check that all symbols are in the GENCODE data
  drivers <- target.gene.prediction.package::read_tibble(driversFile)$V1
  # check that all symbols are in the GENCODE data
  drivers_NotFound <- drivers[drivers %ni% target.gene.prediction.package::TSSs$symbol]
  if(length(drivers_NotFound) > 0){
    message(dplyr::n_distinct(drivers_NotFound),
            " provided driver gene(s) from '", basename(driversFile), "' are not found in the GENCODE gene symbols and therefore cannot be considered.",
            "\nUnknown genes: ", paste(drivers_NotFound, sep="", collapse=", "))
  }

  # ======================================================================================================
  # #### 1) CELL TYPE ENRICHMENT ####
  cat("1) Cell type enrichment...\n")

  specificity_enriched_celltypes <- DHSs %>%
    # Fisher enrichment test of variants in upper-quartile cell-type-specificic H3K27ac marks in DHSs
    dplyr::filter(Method == "specificity",
                  Mark == "H3K27ac",
                  decile == 10) %>%
    target.gene.prediction.package::bed_fisher_grouped(
      bedA = .,
      bedA_groups = "CellType",
      bedB = variants,
      genome = target.gene.prediction.package::ChrSizes,
      # filter for effect and significance
      estimate > 2,
      p.value < 0.05
    ) %>%
    # Extract enriched cell types
    dplyr::pull(CellType) %>%
    {dplyr::filter(target.gene.prediction.package::DHSs_metadata, code %in% .)}

  cat("Enriched cell types: ", specificity_enriched_celltypes$code)

  # Filter to annotations for enriched tissue(s)
  contact_enriched <- target.gene.prediction.package::contact_metadata %>%
    dplyr::filter(Tissue %in% specificity_enriched_celltypes$tissue)
  DHSs_enriched <- target.gene.prediction.package::DHSs_metadata %>%
    dplyr::filter(tissue %in% specificity_enriched_celltypes$tissue)

  # ======================================================================================================
  # #### 2a) GENE-LEVEL INPUTS ####
  # (Intersection with list of annotations, expression profiles, etc...)
  cat("2a) Annotating genes...\n")

  # intersect with list of genomic annotations
  gene_annotations <- target.gene.prediction.package::TSSs %>%
    target.gene.prediction.package::intersect_DHSs(enst)

  # bind and widen all gene-level annotations
  weighted_gene_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = "enst",
    annotation.level = "g",
    gene_annotations
  )

  # ======================================================================================================
  # #### 2b) VARIANT-LEVEL INPUTS ####
  # (Intersection with list of annotations, enhancers, GWAS statistics(?), etc...)
  cat("2b) Annotating", trait, "variants...\n")

  # intersect with list of genomic annotations
  variant_annotations <- variants %>%
    target.gene.prediction.package::intersect_DHSs(variant)

  # calculate n genes near each variant
  variant_n_genes <- variants %>%
    # Get genes within xMb of each variant
    valr::bed_slop(both = variant_to_gene_max_distance, genome = target.gene.prediction.package::ChrSizes, trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs, suffix = c("", ".TSS")) %>%
    dplyr::count(variant) %>%
    dplyr::transmute(variant,
                     annotation.name = paste0("n_genes_within_", variant_to_gene_max_distance/1e6, "Mb_of_variant"),
                     annotation.value = n,
                     annotation.weight = 1)

  # bind and widen all variant-level annotations
  weighted_variant_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = "variant",
    annotation.level = "v",
    variant_annotations,
    variant_n_genes
  )

  # ======================================================================================================
  # #### 2c) GENE-X-VARIANT-LEVEL INPUTS ####
  # (HiC interaction, distance, etc...)
  # The package will consider all genes as potential targets of a variants in CS's whose ranges are within 2Mb of
  cat("2c) Annotating gene x", trait, "variant pairs...\n")

  # TSS distance scoring method (these pairings are the only ones to consider)
  gene_x_variant_distance <- variants %>%
    # Get genes within variant_to_gene_max_distance of each variant
    valr::bed_slop(both = variant_to_gene_max_distance,
                   genome = target.gene.prediction.package::ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs,
                        suffix = c(".variant", ".TSS")) %>%
    dplyr::group_by(variant.variant) %>%
    # Calculate inverse of the absolute bp distance for each variant-gene pair
    dplyr::mutate(pair.distance = abs((start.variant + variant_to_gene_max_distance) - start.TSS),
                  pair.inverse_distance = 1/pair.distance,
                  # ranking transcript TSSs (if two transcript TSSs are equidistant to the variant, they will receive the same, lower rank)
                  pair.distance_rank = rank(pair.distance, ties.method = "min")) %>%
    dplyr::select(variant = variant.variant,
                  cs = cs.variant,
                  enst = enst.TSS,
                  pair.inverse_distance,
                  pair.distance_rank)

  gene_x_variant_distance_score <- gene_x_variant_distance %>%
    dplyr::transmute(variant, cs, enst,
                     annotation.name = "inverse_distance_within_max_distance",
                     annotation.value = pair.inverse_distance,
                     annotation.weight = 1)

  # variant-TSS distance rank method
  gene_x_variant_distance_rank <- gene_x_variant_distance %>%
    dplyr::transmute(variant, cs, enst,
                     annotation.name = "distance_rank_within_max_distance",
                     annotation.value = pair.distance_rank,
                     annotation.weight = 1)

  # intersect loop ends, by cell type, with user-provided variants and gene TSSs
  # (finds interaction loops with a variant at one end and a TSS at the other)
  gene_x_variant_contact <- contact %>%
    purrr::map(~ intersect_BEDPE(
      # ! For mutually exclusive intersection with HiC ranges, make variant intervals 1bp long, equal to the end position
      SNPend = variants %>% dplyr::mutate(start = end),
      TSSend = target.gene.prediction.package::TSSs,
      bedpe = .) %>%
        dplyr::transmute(variant,
                         cs,
                         enst,
                         InteractionID,
                         annotation.value = score,
                         annotation.weight = 1)) %>%
    dplyr::bind_rows(.id = "annotation.name") %>%
    # Make sure all interactions are within 2Mb - hard filter, ignore everything higher
    dplyr::inner_join(gene_x_variant_distance %>% dplyr::select(variant, cs, enst))
  ## ~10% of variant-TSS interactions indicated by the contact data
  ## are further than 2Mb apart and are thus eliminated

  # nearest gene TSS method
  gene_x_variant_closest <- gene_x_variant_distance_rank %>%
    dplyr::filter(annotation.value == 1) %>%
    dplyr::transmute(variant,
                     cs,
                     enst = enst,
                     annotation.name = "nearest",
                     annotation.value = 1,
                     annotation.weight = 1)

  # bind and widen all gene-x-variant-level annotations
  weighted_gene_x_variant_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = c("variant", "cs", "enst"),
    annotation.level = "gxv",
    gene_x_variant_distance_score,
    gene_x_variant_contact,
    gene_x_variant_closest
  )

  # ======================================================================================================
  # #### 2d) GENE-X-CS-LEVEL INPUTS ####
  gene_x_cs_contact <- gene_x_variant_contact %>%
    dplyr::group_by(cs, enst) %>%
    dplyr::mutate(n_gene_x_cs_interactions = dplyr::n_distinct(InteractionID),
                  InteractionID = paste0(unique(InteractionID), collapse = ",")) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(cs,
                     enst,
                     InteractionID,
                     annotation.name = "multi_contact",
                     annotation.value = n_gene_x_cs_interactions,
                     annotation.weight = 1) %>%
    dplyr::distinct()

  weighted_gene_x_cs_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = c("cs", "enst"),
    annotation.level = "gxc",
    gene_x_cs_contact
  )

  # ======================================================================================================
  # #### 3) ALL INPUTS ####
  # Master variant-gene pair table
  # -> wide-format (one row per gene-variant pair, one column per annotation)
  # -> only variant-gene combinations with at least one pair annotation (HiChIP interaction, nearest or within 2Mb) are included
  # -> pair ID columns: | variant | enst |
  # -> annotation columns: | g_* | v_* | gxv_* | gxc_*
  cat("3) Generating master table of gene x", trait, "variant pairs, with all annotation levels (genes, variants, gene-variant pairs, gene-credible set pairs)...\n")
  master <- variants %>%
    dplyr::select(variant, cs) %>%
    # Get all genes with TSSs within variant_to_gene_max_distance of every variant
    dplyr::right_join(  weighted_gene_x_variant_annotations ) %>%
    dplyr::left_join(   weighted_gene_x_cs_annotations      ) %>%
    dplyr::left_join(   weighted_variant_annotations        ) %>%
    dplyr::left_join(   weighted_gene_annotations           ) %>%
    # Replace all NAs in numeric annotation value columns with 0
    dplyr::mutate_if(is.numeric, tidyr::replace_na, replace = 0)

  # ======================================================================================================
  # #### 4) SCORING ####
  # -> Score enriched feature annotations
  predictions_full <- master %>%
    # Add drivers column
    dplyr::left_join(target.gene.prediction.package::TSSs %>% dplyr::select(enst, symbol)) %>%
    dplyr::mutate(driver = symbol %in% drivers) %>%
    # Filtering to only consider CSs with at least 1 driver gene within variant_to_gene_max_distance
    dplyr::group_by(cs) %>%
    dplyr::group_modify(~ .x %>%
                          dplyr::mutate(cs_n_drivers = .x %>%
                                          dplyr::filter(driver) %>%
                                          dplyr::pull(symbol) %>%
                                          dplyr::n_distinct())) %>%
    dplyr::ungroup() %>%
    dplyr::filter(cs_n_drivers > 0)

  predictions <- predictions_full %>%
    # get score, prediction (T/F), and max_score per CS
    dplyr::group_by(cs) %>%
    dplyr::mutate(score = rowSums(dplyr::across(where(is.numeric))),
                  prediction = score > 150,
                  max_score = score == max(score)) %>%
    dplyr::ungroup() %>%
    dplyr::select(cs, variant, enst, symbol, score, prediction, max_score, driver)

  # save output predictions table
  predictions %>%
    dplyr::select(cs, variant, symbol, score) %>%
    write.table(outPredictions,
                quote = F, row.names = F, sep = "\t")

  # ======================================================================================================
  # #### 5) PRECISION-RECALL ####
  # Generate PR curves (model performance metric)
  PR <- target.gene.prediction.package::get_PR(predictions,
                                               score, max_score, prediction)

  # PR of each individual annotation (columns of master)
  PR_all <- target.gene.prediction.package::get_PR(predictions_full,
                                                   gxv_inverse_distance_within_max_distance:gene_DHSs_specificity_H3K27ac_E129_bin10) %>%
    # add annotation level info
    dplyr::mutate(level = sub("_.*", "", prediction_type))


  pdf(outPR, height = 10, onefile = T)
  PR %>%
  PR_all %>%
    dplyr::select(prediction_type, AUC, level) %>%
    dplyr::distinct() %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(prediction_type, AUC),
                        y = AUC,
                        fill = level)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(x = "Predictor", y = "PR AUC")
  PR_all %>%
    dplyr::filter(.threshold %ni% c(0, Inf, -Inf) & recall != 1,
                  grepl("gxv_", prediction_type)) %>%
    dplyr::mutate(line = dplyr::n() > 1,
                  point = dplyr::n() == 1) %>%
    ggplot2::ggplot(ggplot2::aes(x = precision,
                                 y = recall,
                                 colour = prediction_type)) +
    ggplot2::geom_line(data = . %>% dplyr::filter(line)) +
    ggplot2::geom_point(data = . %>% dplyr::filter(point)) +
    ggplot2::theme(legend.position = "none")
  dev.off()

  # ======================================================================================================
  # #### 5) XGBoost model training? ####
  # XGBoost input
}
