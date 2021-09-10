#' Predict target genes of a list of fine-mapped non-coding variants for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param varfile A header-less BED file of fine-mapped trait-relevant variants, with metadata columns for the variant's name and the credible set it belongs to (chrom, start, stop, variant, cs)
#' @param trait The name of the trait of interest.
#' @param tissue The tissue(s) of action for the trait.
#' @param outDir The output directory in which to save the predictions. Default is current directory.
#' @param contactDir The directory containing contact data.
#' @param variant_to_gene_max_distance The maximum absolute distance (bp) across which variant-gene pairs are considered. Default is 2Mb.
#'
#' @return A file of variant-gene pair predictions, with associated scores, saved in the given output directory.
#' @export
predict_target_genes <- function(varfile,
                                 trait = NULL,
                                 tissue = NULL,
                                 outDir = ".",
                                 contactDir = "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/3C/",
                                 variant_to_gene_max_distance = 2e6){

  # silence "no visible binding" NOTE for data variables
  . <- variant <- enst <- start.variant <- start.TSS <- variant.variant <- enst.TSS <- score <- annotation <- estimate <-
    p.value <- enst.query <- annotation.annotation <- variant.query <- n <- end <- pair.score <-
    NULL

  # define the outfile
  outfile <- paste0(outDir,"/") %>% { if(!is.null(trait)) paste0(., trait, "_") else . } %>% paste0("target_gene_predictions.tsv")

  # import the variants
  variants <- target.gene.prediction.package::import_BED(varfile,
                                                         metadata_cols = c("variant", "cs"))

  # import the contact data
  contact <- load_contact_data(contactDir = contactDir)

  # for testing:
  # variants=BCVars %>% dplyr::select(chrom:cs); trait="BC"; outDir="."; contactDir="/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/3C/"; variant_to_gene_max_distance=2e6

  # list annotations to be used by intersect_annotations() (must be a list of BED tibbles)
  annotations <- list(
      DHSs = target.gene.prediction.package::DHSs %>%
        target.gene.prediction.package::recursively_bind_rows(nest_names = c("Method", "Mark", "CellType", "Bin")) %>%
        dplyr::mutate(annotation.description = paste(Bin, "of DHSs binned by", Method, "of", Mark, "annotations in",
                                                     target.gene.prediction.package::DHSs_metadata$name[target.gene.prediction.package::DHSs_metadata==CellType]))
  )
  ## TODO: fix the annotations.rda problem (save online somewhere / GitHub LFS?)
  ## TODO: change `annotations` back to `target.gene.prediction.package::annotations` in this script

  # ======================================================================================================
  # #### 1) CELL TYPE ENRICHMENT ####
  cat("1) Cell type enrichment...\n")

  specificity_enriched_celltypes <- target.gene.prediction.package::DHSs %>%
    target.gene.prediction.package::recursively_bind_rows(nest_names = c("Method", "Mark", "CellType", "Bin")) %>%
    # Fisher enrichment test of variants in upper-quartile cell-type-specificic H3K27ac marks in DHSs
    dplyr::filter(Method == "specificity",
                  Mark == "H3K27ac",
                  Bin == "bin10") %>%
    target.gene.prediction.package::bed_fisher_grouped(
      bedA = .,
      bedA_groups = "CellType",
      bedB = variants,
      genome = target.gene.prediction.package::ChrSizes,
      # filter for effect and significance
      estimate > 2,
      p.value < 0.05
    ) %>%
    # Extract enriched cell types - specificity | H3K27ac | {{ CellType }} | quartiles | 4
    dplyr::pull(CellType) %>%
    {dplyr::filter(target.gene.prediction.package::DHSs_metadata, code %in% .)}

  cat("Enriched cell types: ", specificity_enriched_celltypes$code)

  # ======================================================================================================
  # #### 2a) GENE-LEVEL INPUTS ####
  # (Intersection with list of annotations, expression profiles, etc...)
  cat("2a) Annotating genes...\n")

  # intersect with list of genomic annotations
  gene_annotations <- target.gene.prediction.package::TSSs %>%
    target.gene.prediction.package::intersect_annotations() %>%
    dplyr::transmute(enst = enst.query,
                     annotation.name = annotation.annotation,
                     annotation.description = annotation.description.annotation,
                     annotation.value = 1,
                     annotation.weight = 1)

  # bind and widen all gene-level annotations
  weighted_gene_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = "enst",
    annotation_level = "gene",
    gene_annotations
  )

  # ======================================================================================================
  # #### 2b) VARIANT-LEVEL INPUTS ####
  # (Intersection with list of annotations, enhancers, GWAS statistics(?), etc...)
  cat("2b) Annotating", trait, "variants...\n")

  # intersect with list of genomic annotations
  variant_annotations <- variants %>%
    target.gene.prediction.package::intersect_annotations() %>%
    dplyr::transmute(variant = variant.query,
                     annotation.name = annotation.annotation,
                     annotation.value = 1,
                     annotation.weight = 1)

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
    annotation_level = "variant",
    variant_annotations,
    variant_n_genes
  )

  # ======================================================================================================
  # #### 2c) GENE-X-VARIANT-LEVEL INPUTS ####
  # (HiC interaction, distance, etc...)
  # The package will consider all genes as potential targets of a variants in CS's whose ranges are within 2Mb of
  cat("2c) Annotating gene x", trait, "variant pairs...\n")

  # intersect HiC ends, by cell type, with user-provided variants and gene TSSs
  # (finds interaction loops with a variant at one end and a TSS at the other)
  gene_x_variant_contact <- contact %>%
    purrr::map(~ intersect_BEDPE(
      # ! For the purpose of mutually exclusive intersection with HiC ranges, make variant intervals 1bp long, equal to the end position
      SNPend = variants %>% dplyr::mutate(start = end),
      TSSend = target.gene.prediction.package::TSSs,
      bedpe = .) %>%
        dplyr::transmute(variant,
                         cs,
                         enst,
                         InteractionID,
                         annotation.value = 1,
                         annotation.weight = 1)) %>%
    dplyr::bind_rows(.id = "annotation.name")

  # nearest gene TSS method
  gene_x_variant_closest <- valr::bed_closest(x = variants,
                                              y = target.gene.prediction.package::TSSs,
                                              suffix = c("", ".TSS")) %>%
    dplyr::distinct() %>%
    dplyr::transmute(variant,
                     enst = enst.TSS,
                     annotation.name = "nearest",
                     annotation.value = 1,
                     annotation.weight = 1)

  # TSS distance scoring method
  gene_x_variant_distance_score <- variants %>%
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
                     annotation.value = pair.score,
                     annotation.weight = 1)

  # bind and widen all gene-x-variant-level annotations
  weighted_gene_x_variant_annotations <- target.gene.prediction.package::bind_and_weight_and_widen_annotations(
    id_cols = c("variant", "enst", "InteractionID"),
    annotation_level = "gxv",
    gene_x_variant_contact,
    gene_x_variant_closest,
    gene_x_variant_distance_score
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
    id_cols = c("cs", "enst", "InteractionID"),
    annotation_level = "gxc",
    gene_x_cs_contact
  )

  # ======================================================================================================
  # #### 3) ALL INPUTS ####
  # Master variant-gene pair table
  # -> wide-format (one row per gene-variant pair, one column per annotation)
  # -> only variant-gene combinations with at least one pair annotation (HiChIP interaction, nearest or within 2Mb) are included
  # -> pair ID columns: | variant | enst |
  # -> annotation columns: | pair_* | gene_* | variant_* |
  cat("3) Generating master table of gene x", trait, "variant pairs, with all annotation levels (genes, variants, gene-variant pairs)...\n")
  master <- variants %>%
    dplyr::select(variant, cs) %>%
    dplyr::right_join(  weighted_gene_x_variant_annotations ) %>%
    dplyr::left_join(   weighted_gene_x_cs_annotations      ) %>%
    dplyr::left_join(   weighted_variant_annotations        ) %>%
    dplyr::left_join(   weighted_gene_annotations           )
  return(master)

  # ======================================================================================================
  # #### 4) SCORING ####
  # -> Score enriched feature annotations

  # ======================================================================================================
  # #### 5) PRECISION-RECALL ####
  # Generate PR curves (model performance metric)

  # predictive value of each individual feature



  # ======================================================================================================
  # #### 5) XGBoost model training? ####

}
