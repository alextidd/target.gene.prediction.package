#' Predict target genes of a list of fine-mapped non-coding variants for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param varfile A header-less BED file of fine-mapped trait-relevant variants, with metadata columns for the variant's name and the credible set it belongs to (chrom, start, stop, variant, cs)
#' @param trait The name of the trait of interest.
#' @param tissue The tissue(s) of action for the trait.
#' @param outdir The output directory in which to save the predictions. Default is current directory.
#' @param variant_to_gene_max_distance The maximum absolute distance (bp) across which variant-gene pairs are considered. Default is 2Mb.
#'
#' @return A file of variant-gene pair predictions, with associated scores, saved in the given output directory.
#' @export
predict_target_genes <- function(varfile, trait = NULL, tissue = NULL, outdir = ".", variant_to_gene_max_distance = 2e6){

  # silence "no visible binding" NOTE for data variables
  . <- variant <- enst <- start.variant <- start.TSS <- variant.variant <- enst.TSS <- score <- annotation <- estimate <-
    p.value <- enst.query <- annotation.annotation <- variant.query <- n <- end <- pair.score <-
    NULL

  # define the outfile
  outfile <- paste0(outdir,"/") %>% { if(!is.null(trait)) paste0(., trait, "_") else . } %>% paste0("target_gene_predictions.tsv")

  # import the variants
  variants <- target.gene.prediction.package::import_BED(varfile,
                                                         metadata_cols = c("variant", "cs"))

  # ======================================================================================================
  # #### 1) CELL TYPE ENRICHMENT ####
  cat("1) Cell type enrichment...\n")

  specificity_enriched_annotations <- target.gene.prediction.package::annotations[["DHSs"]] %>%
    # Fisher enrichment test of variants in upper-quartile cell-type-specificic H3K27ac marks in DHSs
    dplyr::filter(grepl("DHSs_specificity_H3K27ac_.*_quartiles_4", annotation)) %>%
    target.gene.prediction.package::bed_fisher_grouped(
      bedA = .,
      bedA_groups = "annotation",
      bedB = variants,
      genome = target.gene.prediction.package::ChrSizes,
      estimate > 1.5,
      p.value < 0.05
    ) %>%
    # Extract enriched cell types from annotation names - DHSs_specificity_H3K27ac_{{ CellType }}_quartiles_4
    dplyr::mutate(code = annotation %>%
                    sub("DHSs_specificity_H3K27ac_", "", .) %>%
                    sub("_quartiles_4", "", .)) %>%
    dplyr::left_join(target.gene.prediction.package::annotations_metadata[["DHSs"]])

  target.gene.prediction.package::annotations[["TFBSs"]] %>%
    # Annotations in the tissues of interest
    dplyr::mutate(code = annotation %>%
                    sub(".*_", "", .)) %>%
    dplyr::left_join(target.gene.prediction.package::annotations_metadata[["TFBSs"]]) %>%
    dplyr::filter(tissue %in% specificity_enriched_annotations$tissue) %>%
    # Fisher enrichment test of variants in TFBSs of tissue type(s) of interest + filter for effect and significance
    target.gene.prediction.package::bed_fisher_grouped(
      bedA = .,
      bedA_groups = "annotation",
      bedB = variants,
      genome = target.gene.prediction.package::ChrSizes,
      estimate > 1.5,
      p.value < 0.05
    ) %>%
    # Extract enriched cell types from annotation names - DHSs_specificity_H3K27ac_{{ CellType }}_quartiles_4
    dplyr::mutate(code = annotation %>%
                    sub("DHSs_specificity_H3K27ac_", "", .) %>%
                    sub("_quartiles_4", "", .)) %>%
    dplyr::left_join(target.gene.prediction.package::annotations_metadata[["DHSs"]])


  # ======================================================================================================
  # #### 2a) GENE-LEVEL INPUTS ####
  # (Intersection with list of annotations, expression profiles, etc...)
  cat("2a) Annotating genes...\n")

  # intersect with list of genomic annotations
  gene_annotations <- target.gene.prediction.package::TSSs %>%
    target.gene.prediction.package::intersect_annotations_list() %>%
    dplyr::transmute(enst = enst.query,
                     annotation.name = annotation.annotation,
                     annotation.value = "TRUE")

  # bind and widen all gene-level annotations
  gene_annotations <- target.gene.prediction.package::bind_and_widen_annotations(
    id_cols = "enst",
    annotation.level = "gene",
    gene_annotations
    )

  # ======================================================================================================
  # #### 2b) VARIANT-LEVEL INPUTS ####
  # (Intersection with list of annotations, enhancers, GWAS statistics(?), etc...)
  cat("2b) Annotating", trait, "variants...\n")

  # intersect with list of genomic annotations
  variant_annotations <- variants %>%
    target.gene.prediction.package::intersect_annotations_list() %>%
    dplyr::transmute(variant = variant.query,
                     annotation.name = annotation.annotation,
                     annotation.value = "TRUE")

  # calculate n genes near each variant
  variant_n_genes <- variants %>%
    # Get genes within 2Mb of each variant
    valr::bed_slop(both = variant_to_gene_max_distance, genome = target.gene.prediction.package::ChrSizes, trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs, suffix = c("", ".TSS")) %>%
    dplyr::count(variant) %>%
    dplyr::transmute(variant,
                     annotation.name = paste0("n_genes_within_", variant_to_gene_max_distance/1e6, "Mb_of_variant"),
                     annotation.value = as.character(n))

  # bind and widen all variant-level annotations
  variant_annotations <- target.gene.prediction.package::bind_and_widen_annotations(
    id_cols = "variant",
    annotation.level = "variant",
    variant_annotations,
    variant_n_genes
  )

  # ======================================================================================================
  # #### 2c) PAIR-LEVEL INPUTS ####
  # (HiChIP interaction, distance, etc...)
  cat("2c) Annotating gene x", trait, "variant pairs...\n")

  # intersect HiChIP ends with user-provided variants and gene TSSs (finds interaction loops with a variant at one end and a TSS at the other)
  pair_hichip <- variants %>%
    # ! For the purpose of mutually exclusive intersection with HiChIP ranges, make variant intervals 1bp long, equal to the end position
    dplyr::mutate(start = end) %>%
    target.gene.prediction.package::intersect_BEDPE(SNPend = .,
                                                    TSSend = target.gene.prediction.package::TSSs,
                                                    bedpe = target.gene.prediction.package::hichip) %>%
    dplyr::transmute(variant,
                     enst,
                     annotation.name = "HiChIP",
                     annotation.value = "TRUE")

  # nearest gene TSS method
  pair_closest <- valr::bed_closest(x = variants,
                                    y = target.gene.prediction.package::TSSs,
                                    suffix = c("", ".TSS")) %>%
    dplyr::distinct() %>%
    dplyr::transmute(variant,
                     enst = enst.TSS,
                     annotation.name = "nearest",
                     annotation.value = "TRUE")

  # TSS distance scoring method
  pair_distance_score <- variants %>%
    # Get genes within 1Mb of each variant
    valr::bed_slop(both = variant_to_gene_max_distance,
                   genome = target.gene.prediction.package::ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs,
                        suffix = c(".variant", ".TSS")) %>%
    # Calculate inverse of the absolute bp distance for each variant-gene pair
    dplyr::mutate(pair.score = 1/(abs((start.variant + variant_to_gene_max_distance) - start.TSS))) %>%
    dplyr::transmute(variant = variant.variant,
                     enst = enst.TSS,
                     annotation.name = "inverse_distance_within_2Mb",
                     annotation.value = as.character(pair.score))

  # bind and widen all pair-level annotations
  pair_annotations <- target.gene.prediction.package::bind_and_widen_annotations(
    id_cols = c("variant", "enst"),
    annotation.level = "pair",
    pair_hichip,
    pair_closest,
    pair_distance_score
  )

  # ======================================================================================================
  # #### 3) ALL INPUTS ####
  # Master variant-gene pair table
  # -> wide-format (one row per pair, one column per annotation)
  # -> only variant-gene combinations with at least one pair annotation (HiChIP interaction, nearest or within 2Mb) are included
  # -> pair ID columns: | variant | enst |
  # -> annotation columns: | pair_* | gene_* | variant_* |
  cat("3) Generating master table of gene x", trait, "variant pairs, with all annotation levels (genes, variants, gene-variant pairs)...\n")
  master <- pair_annotations %>%
    dplyr::left_join(variant_annotations) %>%
    dplyr::left_join(gene_annotations)
  return(master)

  # ======================================================================================================
  # #### 4) SCORING ####
  # -> Score enriched feature annotations

}
