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
  outBase <- paste0(outDir,"/") %>% { if(!is.null(trait)) paste0(., trait, "_") else . }
  outAnnotations <- paste0(outBase, "target_gene_annotations.tsv")
  outPredictions <- paste0(outBase, "target_gene_predictions.tsv")
  outPR <- paste0(outBase, "PrecisionRecall.pdf")

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
  weight <- dplyr::tibble(
    v_n_genes = 1,
    gxv_distance_score = 1,
    gxv_distance_rank = 1,
    gxv_closest = 1,
    gxv_contact = 2,
    gxc_multicontact = 3
  )

  # ======================================================================================================
  # #### 1) CELL TYPE ENRICHMENT ####
  cat("1) Cell type enrichment...\n")
  enriched_celltypes <- DHSs %>%
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
    {dplyr::filter(DHSs_metadata, code %in% .)}

  enriched_tissues <- DHSs_metadata %>%
    dplyr::filter(tissue %in% enriched_celltypes$tissue)

  cat("Enriched cell types: ", enriched_celltypes$code, "\n")
  cat("Enriched tissues:", unique(enriched_celltypes$tissue), "\n")

  # Subset DHS annotations
  enriched_DHSs <- DHSs %>%
    dplyr::filter(CellType %in% enriched_tissues$code)
  enriched_contact_elements <- contact_metadata %>%
    dplyr::filter(Tissue %in% enriched_tissues$tissue) %>%
    dplyr::pull(list_element)

  # ======================================================================================================
  # #### 2) ENHANCER VARIANTS ####
  cat("1) Finding enhancer variants...\n")
  enh_variants <- variants %>%
    target.gene.prediction.package::bed_intersect_left(DHSs, keepBcoords = F, keepBmetadata = F)

  # ======================================================================================================
  # #### 3a) ENHANCER-LEVEL INPUTS ####
  # (Intersection with list of annotations, enhancers, GWAS statistics(?), etc...)
  cat("2b) Annotating", trait, "enhancer variants...\n")

  # intersect with list of genomic annotations
  v_annotations <- enh_variants %>%
    target.gene.prediction.package::intersect_DHSs(
      variant,
      annots = enriched_DHSs) %>% # %>% dplyr::filter(Method == "signal")) %>%
    dplyr::mutate(annotation.weight = 1)

  # calculate n genes near each variant
  v_n_genes <- enh_variants %>%
    # Get genes within xMb of each variant
    valr::bed_slop(both = variant_to_gene_max_distance, genome = target.gene.prediction.package::ChrSizes, trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs, suffix = c("", ".TSS")) %>%
    dplyr::count(variant) %>%
    dplyr::transmute(variant,
                     annotation.name = paste0("inverse_n_genes_within_", variant_to_gene_max_distance/1e6, "Mb_of_variant"),
                     annotation.value = 1/n,
                     annotation.weight = weight$v_n_genes)

  # bind and widen all variant-level annotations
  weighted_v_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = "variant",
    annotation.level = "v",
    v_annotations,
    v_n_genes
  )

  # ======================================================================================================
  # #### 3b) GENE-LEVEL INPUTS ####
  # (Intersection with list of annotations, expression profiles, etc...)
  cat("2a) Annotating genes...\n")

  # intersect with DHSs
  g_annotations <- target.gene.prediction.package::TSSs %>%
    target.gene.prediction.package::intersect_DHSs(
      enst,
      annots = enriched_DHSs) %>% # %>% dplyr::filter(Method == "signal")) %>%
    dplyr::mutate(annotation.weight = 1)

  # bind and widen all gene-level annotations
  weighted_g_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = "enst",
    annotation.level = "g",
    g_annotations
  )

  # ======================================================================================================
  # #### 3c) GENE-X-VARIANT-LEVEL INPUTS ####
  # (HiC interaction, distance, etc...)
  # The package will consider all genes as potential targets of a variant in CS's whose ranges are within 2Mb of
  cat("2c) Annotating gene x", trait, "variant pairs...\n")

  # TSS distance scoring method (these pairings are the only ones to consider)
  gxv_distance <- enh_variants %>%
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
                  pair.inverse_distance_rank = 1/rank(pair.distance, ties.method = "min")) %>%
    dplyr::select(variant = variant.variant,
                  cs = cs.variant,
                  enst = enst.TSS,
                  pair.inverse_distance,
                  pair.inverse_distance_rank)

  # variant-TSS inverse distance score method
  gxv_distance_score <- gxv_distance %>%
    dplyr::transmute(variant, cs, enst,
                     annotation.name = "inverse_distance_within_max_distance",
                     annotation.value = pair.inverse_distance,
                     annotation.weight = weight$gxv_distance_score)

  # variant-TSS distance rank method
  gxv_distance_rank <- gxv_distance %>%
    dplyr::transmute(variant, cs, enst,
                     annotation.name = "inverse_distance_rank_within_max_distance",
                     annotation.value = pair.inverse_distance_rank,
                     annotation.weight = weight$gxv_distance_rank)

  # closest variant-TSS method
  gxv_closest <- gxv_distance %>%
    dplyr::filter(pair.inverse_distance_rank == 1) %>%
    dplyr::transmute(variant, cs, enst,
                     annotation.name = "closest",
                     annotation.value = 1,
                     annotation.weight = weight$gxv_closest)

  # intersect loop ends, by cell type, with enhancer variants and gene TSSs
  # (finds interaction loops with a variant at one end and a TSS at the other)
  gxv_contact <- contact %>%
    purrr::map(~ intersect_BEDPE(
      # ! For mutually exclusive intersection with HiC ranges, make variant intervals 1bp long, equal to the end position
      SNPend = enh_variants %>% dplyr::mutate(start = end),
      TSSend = target.gene.prediction.package::TSSs,
      bedpe = .) %>%
        dplyr::transmute(variant,
                         cs,
                         enst,
                         InteractionID,
                         annotation.value = score)) %>%
    dplyr::bind_rows(.id = "annotation.name") %>%
    dplyr::mutate(
      annotation.weight = dplyr::case_when(annotation.name %in% enriched_contact_elements ~ 2 * weight$gxv_contact,
                                           TRUE ~ weight$gxv_contact)) %>%
    # Make sure all interactions are within 2Mb - hard filter, ignore everything further
    dplyr::inner_join(gxv_distance %>% dplyr::select(variant, cs, enst))
  ## ~10% of variant-TSS interactions indicated by the contact data
  ## are further than 2Mb apart and are thus eliminated

  # bind and widen all gene-x-variant-level annotations
  weighted_gxv_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = c("variant", "cs", "enst"),
    annotation.level = "gxv",
    gxv_closest,
    gxv_distance_score,
    gxv_distance_rank,
    gxv_contact
  )

  # ======================================================================================================
  # #### 3d) GENE-X-CS-LEVEL INPUTS ####
  gxc_multicontact <- gxv_contact %>%
    # multicontact statistics within each variant-x-cs-x-experiment combination
    dplyr::group_by(cs, enst, annotation.name) %>%
    dplyr::mutate(
      # number of gxc loops within each experiment
      n_gxc_contacts = dplyr::n_distinct(InteractionID),
      # sum of the scores of gxc loops within each experiment
      sum_gxc_contacts = sum(annotation.value),
      # collapse Interaction IDs
      InteractionID = paste0(unique(InteractionID), collapse = ",")
    ) %>%
    # This was taken out (only one instance in HaCaT.stimulated_CHiC, negligible improvement)
    # # multicontact statistics within each variance-x-cs-x-celltype combination (if there is >1 assay per celltype, combine and assign max scores per interaction)
    # dplyr::group_by(cs, enst, CellType) %>%
    # dplyr::mutate(
    #   max_sum_gxc_contacts = max(sum_gxc_contacts)
    # ) %>%
    dplyr::ungroup()  %>%
    dplyr::transmute(cs,
                     enst,
                     InteractionID,
                     annotation.name = paste0("multicontact_", annotation.name),
                     annotation.value = sum_gxc_contacts,
                     # weighting (more for enriched tissues)
                     annotation.weight = dplyr::case_when(annotation.name %in% enriched_contact_elements ~ 2 * weight$gxc_multicontact,
                                                          TRUE ~ weight$gxc_multicontact)) %>%
    dplyr::distinct()

  weighted_gxc_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = c("cs", "enst"),
    annotation.level = "gxc",
    gxc_multicontact
  )

  # ======================================================================================================
  # #### 4) ALL INPUTS ####
  # Master variant-gene pair table
  # -> wide-format (one row per gene-variant pair, one column per annotation)
  # -> only variant-gene combinations with at least one pair annotation (HiChIP interaction, nearest or within 2Mb) are included
  # -> pair ID columns: | variant | enst |
  # -> annotation columns: | g_* | v_* | gxv_* | gxc_*
  cat("3) Generating master table of gene x", trait, "variant pairs, with all annotation levels (genes, variants, gene-variant pairs, gene-credible set pairs)...\n")
  master <- enh_variants %>%
    dplyr::select(variant, cs) %>%
    # Get annotations for all genes with TSSs within variant_to_gene_max_distance of every variant
    dplyr::right_join(  weighted_gxv_annotations      ) %>%
    dplyr::left_join(   weighted_gxc_annotations      ) %>%
    dplyr::left_join(   weighted_v_annotations        ) %>%
    dplyr::left_join(   weighted_g_annotations        ) %>%
    # Replace all NAs in numeric annotation value columns with 0
    dplyr::mutate_if(is.numeric, tidyr::replace_na, replace = 0) %>%
    # Add gene symbols and logical driver column (T/F)
    dplyr::left_join(target.gene.prediction.package::TSSs %>% dplyr::select(enst, symbol)) %>%
    dplyr::mutate(driver = symbol %in% drivers)

  predictions <- master %>%
    # Score (row sum)
    dplyr::mutate(score = rowSums(dplyr::across(where(is.numeric)))) %>%
    dplyr::select(cs, symbol, score, driver) %>%
    dplyr::distinct() %>%
    # Get maximum scoring variant per gene-x-CS
    dplyr::group_by(cs, symbol) %>%
    dplyr::filter(score == max(score)) %>%
    # Get maximum scoring gene per CS
    dplyr::group_by(cs) %>%
    dplyr::mutate(max = score == max(score)) %>%
    dplyr::ungroup() %>%
    # Get 'predictions' based on cutoff
    dplyr::mutate(prediction = score > 8)

  # save output tables
  cat("Saving prediction files...")
  master %>%
    write.table(outAnnotations,
                quote = F, row.names = F, sep = "\t")
  predictions %>%
    dplyr::select(cs, variant, symbol, score) %>%
    write.table(outPredictions,
                quote = F, row.names = F, sep = "\t")


  # ======================================================================================================
  # #### 6) PRECISION-RECALL ####
  # Generate PR curves (model performance metric)
  annotation_cols <- names(master)[greplany(c("g_","gxv_","gxc_","v_"),names(master))]

  PR <- predictions %>%
    # Filtering to only consider CSs with at least 1 driver gene within variant_to_gene_max_distance
    dplyr::group_by(cs) %>%
    dplyr::filter(any(driver==T)) %>%
    dplyr::group_by(cs, symbol) %>%
    # PR
    target.gene.prediction.package::get_PR(.,
                                           score, max, prediction)

  # PR of each individual annotation (columns of master)
  PR_all <- target.gene.prediction.package::get_PR(predictions_full,
                                                   annotation_cols) %>%
    # add annotation level info
    dplyr::mutate(level = sub("_.*", "", prediction_type))


  pdf(outPR, height = 10, onefile = T)
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
    dplyr::filter()
    plot_PR() +
    ggplot2::theme(legend.position = "none")
  dev.off()

  # ======================================================================================================
  # #### 7) XGBoost model training? ####
  # XGBoost input
  master %>%
    dplyr::select(label = driver,
                  annotation_cols)


}
