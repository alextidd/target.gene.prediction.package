get_v_level_annotations <- function(open.variants = open_variants,
                                    enriched.DHSs = enriched$DHSs,
                                    variant.to.gene.max.distance = variant_to_gene_max_distance){
  cat("Annotating variants...\n")

  v <- list()

  # intersect with DHS binnings
  v <- v %>%
    intersect_DHSs(open.variants %>% dplyr::select(chrom:variant),
                   enriched.DHSs)

  # calculate n genes near each variant
  v$inv_n_genes <- open.variants %>%
    # Get genes within xMb of each variant
    valr::bed_slop(both = variant.to.gene.max.distance,
                   genome = target.gene.prediction.package::ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs, suffix = c("", ".TSS")) %>%
    dplyr::count(variant) %>%
    dplyr::transmute(variant,
                     value = 1/n)


  # return
  names(v) <- paste0("v_", names(v))
  return(v)
}
