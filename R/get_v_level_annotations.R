get_v_level_annotations <- function(){
  cat("Annotating", trait, "variants...\n")

  v <- list()

  # intersect with list of genomic annotations
  v$annotations <- open_variants %>%
    target.gene.prediction.package::intersect_DHSs(
      query = .,
      variant,
      annots = enriched_DHSs) %>%
    dplyr::mutate(annotation.weight = 1)

  # calculate n genes near each variant
  v$n_genes <- open_variants %>%
    # Get genes within xMb of each variant
    valr::bed_slop(both = variant_to_gene_max_distance, genome = target.gene.prediction.package::ChrSizes, trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs, suffix = c("", ".TSS")) %>%
    dplyr::count(variant) %>%
    dplyr::transmute(variant,
                     annotation.name = paste0("inverse_n_genes_within_", variant_to_gene_max_distance/1e6, "Mb_of_variant"),
                     annotation.value = 1/n,
                     annotation.weight = weight$v$n_genes)

  # return
  return(v)
}
