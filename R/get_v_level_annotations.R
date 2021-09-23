get_v_level_annotations <- function(open.variants = open_variants,
                                    enriched.DHSs = enriched$DHSs,
                                    variant.to.gene.max.distance = variant_to_gene_max_distance,
                                    weight.v = weight$v){
  cat("Annotating variants...\n")

  v <- list()

  # intersect with list of genomic annotations
  v$annotations <- target.gene.prediction.package::intersect_DHSs(
      query = open.variants %>% dplyr::select(-DHS),
      DHSs = enriched.DHSs,
      variant) %>%
    dplyr::mutate(annotation.weight = 1)

  # calculate n genes near each variant
  v$n_genes <- open.variants %>%
    # Get genes within xMb of each variant
    valr::bed_slop(both = variant.to.gene.max.distance,
                   genome = target.gene.prediction.package::ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs, suffix = c("", ".TSS")) %>%
    dplyr::count(variant) %>%
    dplyr::transmute(variant,
                     annotation.name = paste0("inverse_n_genes_within_", variant.to.gene.max.distance/1e6, "Mb_of_variant"),
                     annotation.value = 1/n,
                     annotation.weight = weight.v$n_genes)

  # return
  return(v)
}
