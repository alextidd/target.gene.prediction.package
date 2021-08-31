#' Predict target genes of a list of fine-mapped non-coding SNPs for a trait
#'
#' The master, user-facing function of this package.
#'
#' @param snpfile A header-less BED file of fine-mapped trait-relevant SNPs, with metadata columns for the variant's name and the credible set it belongs to (chrom, start, stop, variant, cs)
#' @param trait The name of the trait of interest
#' @param outdir The output directory in which to save the predictions. If none given, saved in current directory.
#'
#' @return A file of variant-gene pair predictions, with associated scores, saved in the given output directory.
#' @export
predict_target_genes <- function(snpfile, trait = NULL, outdir = "."){

  # silence "no visible binding" NOTE for data variables
  . <- variant <- enst <- start.SNP <- start.TSS <- variant.SNP <- enst.TSS <- score <- NULL

  # Define the outfile
  outfile <- paste0(outdir,"/") %>% { if(!is.null(trait)) paste0(., trait, "_") else . } %>% paste0("target_gene_predictions.tsv")

  # Import the variants
  snps <- target.gene.prediction.package::import_BED(snpfile,
                                                     metadata_cols = c("variant", "cs"))

  # ======================================================================================================
  # GENE-LEVEL INPUTS
  # (Intersection with list of annotations, expression profiles, etc...)
  gene_annotations <- list()
  for(annotation in names(target.gene.prediction.package::annotations)){
    cat("============================================\n",
        "Intersecting genes with", annotation, "\n")
    gene_annotations[[annotation]] <- target.gene.prediction.package::bed_intersect_left(target.gene.prediction.package::TSSs,
                                                                                         target.gene.prediction.package::annotations[[annotation]],
                                                                                         keepBcoords = F) %>%
      dplyr::select(-setdiff(names(target.gene.prediction.package::TSSs), "enst"))
  }

  # ======================================================================================================
  # SNP-LEVEL INPUTS
  # (Intersection with list of annotations, number of genes, GWAS statistics(?), etc...)
  snp_annotations <- list()
  for(annotation in names(target.gene.prediction.package::annotations)){
    cat("============================================\n",
        "Intersecting", trait, "SNPs with", annotation, "\n")
    snp_annotations[[annotation]] <- target.gene.prediction.package::bed_intersect_left(snps,
                                                                                        target.gene.prediction.package::annotations[[annotation]],
                                                                                        keepBcoords = F) %>%
      dplyr::select(-setdiff(names(snps), "variant"))
  }

  n_genes_per_snp <- snps %>%
    # Get genes within 1Mb of each SNP
    valr::bed_slop(both = 1e6, genome = target.gene.prediction.package::ChrSizes, trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs, suffix = c("", ".TSS")) %>%
    dplyr::count(variant)

  # ======================================================================================================
  # PAIR-LEVEL INPUTS
  # (HiChIP interaction, distance, etc...)

  # intersect HiChIP ends with user-provided SNPs and gene TSSs (finds interaction loops with a SNP at one end and a TSS at the other)
  hichip_pairs <- target.gene.prediction.package::intersect_BEDPE(SNPbed = snps,
                                                                  TSSbed = target.gene.prediction.package::TSSs,
                                                                  bedpe = target.gene.prediction.package::hichip) %>%
    dplyr::select(variant, enst)

  # nearest gene TSS method
  closest_pairs <- valr::bed_closest(x = snps,
                                     y = target.gene.prediction.package::TSSs,
                                     suffix = c("", "")) %>%
    dplyr::select(variant, enst)

  # TSS distance scoring method
  distance_scored_pairs <- snps %>%
    # Get genes within 1Mb of each SNP
    valr::bed_slop(both = 1e6, genome = target.gene.prediction.package::ChrSizes, trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs, suffix = c(".SNP", ".TSS")) %>%
    # Calculate inverse of the absolute bp distance for each SNP-gene pair
    dplyr::mutate(score = 1/(abs((start.SNP + 1e6) - start.TSS))) %>%
    dplyr::select(variant = variant.SNP, enst = enst.TSS, score)

}
