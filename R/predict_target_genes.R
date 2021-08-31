#' Predict target genes of a list of fine-mapped non-coding SNPs for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param snpfile A header-less BED file of fine-mapped trait-relevant SNPs, with metadata columns for the variant's name and the credible set it belongs to (chrom, start, stop, variant, cs)
#' @param trait The name of the trait of interest
#' @param tissue The tissue(s) of action for the trait
#' @param outdir The output directory in which to save the predictions. If none given, saved in current directory.
#'
#' @return A file of variant-gene pair predictions, with associated scores, saved in the given output directory.
#' @export
predict_target_genes <- function(snpfile, trait = NULL, tissue = NULL, outdir = "."){

  # silence "no visible binding" NOTE for data variables
  . <- variant <- enst <- start.SNP <- start.TSS <- variant.SNP <- enst.TSS <- score <- NULL

  # define the outfile
  outfile <- paste0(outdir,"/") %>% { if(!is.null(trait)) paste0(., trait, "_") else . } %>% paste0("target_gene_predictions.tsv")

  # import the variants
  snps <- target.gene.prediction.package::import_BED(snpfile,
                                                     metadata_cols = c("variant", "cs"))

  # ======================================================================================================
  # GENE-LEVEL INPUTS
  # (Intersection with list of annotations, expression profiles, etc...)

  # list of genomic annotations
  gene_annotations <- list() ; for(annotation in names(target.gene.prediction.package::annotations)){
    cat("============================================\n",
        "Intersecting gene TSSs with", annotation, "\n")
    gene_annotations[[annotation]] <- target.gene.prediction.package::bed_intersect_left(target.gene.prediction.package::TSSs,
                                                                                         target.gene.prediction.package::annotations[[annotation]],
                                                                                         keepBcoords = F) %>%
      dplyr::select(enst, gene.annotation = annotation)
  } ; gene_annotations <- dplyr::bind_rows(gene_annotations)

  # ======================================================================================================
  # SNP-LEVEL INPUTS
  # (Intersection with list of annotations, enhancers, GWAS statistics(?), etc...)

  # list of genomic annotations
  variant_annotations <- list() ; for(annotation in names(target.gene.prediction.package::annotations)){
    cat("============================================\n",
        "Intersecting", trait, "SNPs with", annotation, "\n")
    variant_annotations[[annotation]] <- target.gene.prediction.package::bed_intersect_left(snps,
                                                                                        target.gene.prediction.package::annotations[[annotation]],
                                                                                        keepBcoords = F) %>%
      dplyr::select(variant, variant.annotation = annotation)
  } ; variant_annotations <- dplyr::bind_rows(variant_annotations)

  # n genes near to the SNP
  variant_n_genes <- snps %>%
    # Get genes within 1Mb of each SNP
    valr::bed_slop(both = 1e6, genome = target.gene.prediction.package::ChrSizes, trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs, suffix = c("", ".TSS")) %>%
    dplyr::count(variant) %>%
    dplyr::transmute(variant,
                     variant.annotation = "n_genes",
                     variant.score = n)

  # ======================================================================================================
  # PAIR-LEVEL INPUTS
  # (HiChIP interaction, distance, etc...)

  # intersect HiChIP ends with user-provided SNPs and gene TSSs (finds interaction loops with a SNP at one end and a TSS at the other)
  pair_hichip <- target.gene.prediction.package::intersect_BEDPE(SNPend = snps,
                                                                  TSSend = target.gene.prediction.package::TSSs,
                                                                  bedpe = target.gene.prediction.package::hichip) %>%
    dplyr::transmute(variant,
                     enst,
                     pair.annotation = "HiChIP")

  # nearest gene TSS method
  pair_closest <- valr::bed_closest(x = snps,
                                    y = target.gene.prediction.package::TSSs,
                                    suffix = c("", ".TSS")) %>%
    dplyr::transmute(variant,
                     enst = enst.TSS,
                     pair.annotation = "nearest")

  # TSS distance scoring method
  pair_distance_scored <- snps %>%
    # Get genes within 1Mb of each SNP
    valr::bed_slop(both = 1e6, genome = target.gene.prediction.package::ChrSizes, trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs, suffix = c(".SNP", ".TSS")) %>%
    # Calculate inverse of the absolute bp distance for each SNP-gene pair
    dplyr::mutate(pair.score = 1/(abs((start.SNP + 1e6) - start.TSS))) %>%
    dplyr::transmute(variant = variant.SNP,
                     enst = enst.TSS,
                     pair.annotation = "inverse_distance_within_1Mb",
                     pair.score)

  # ======================================================================================================
  # ALL INPUTS - Master SNP-gene table
  master <- dplyr::full_join(pair_hichip, pair_closest) %>%
    dplyr::full_join(pair_distance_scored) %>%
    dplyr::left_join(variant_annotations) %>%
    dplyr::left_join(variant_n_genes) %>%
    dplyr::left_join(gene_annotations)

}
