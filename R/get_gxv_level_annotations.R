get_gxv_level_annotations <- function(variants,
                                      txv_master,
                                      enriched){

  gxv <- list()

  # variant-gene closest distance (among all of the gene's transcripts' TSSs)
  gxv$inv_distance_rank <- txv_master %>%
    dplyr::group_by(variant, symbol) %>%
    dplyr::summarise(distance = min(distance)) %>%
    dplyr::group_by(variant) %>%
    dplyr::transmute(
      symbol,
      variant,
      value = 1/rank(distance, ties.method = "min")
    )

  # variants in specific  x genes in specific H3K27ac
  gxv$specific_H3K27ac_closest_specific_genes <- variants %>%
    # intersect with specific_H3K27ac_closest_specific_genes
    bed_intersect_left(enriched$specific_H3K27ac_closest_specific_genes, keepBcoords = F) %>%
    tidyr::pivot_longer(-c(chrom:cs), names_to = "celltype", values_to = "ensg") %>%
    dplyr::filter(!is.na(ensg)) %>%
    dplyr::distinct(variant, celltype, ensg) %>%
    dplyr::mutate(value = 1) %>%
    tidyr::pivot_wider(id_cols = -c(celltype, value),
                       names_from = celltype, values_from = value)

  # REVEL #
  intersect_REVEL <- function(variants, revel_df, ...){
    variants %>%
      dplyr::inner_join(revel_df, by = c("chrom" = "chrom", "end" = "position")) %>%
      tidyr::separate_rows(ensgs) %>%
      dplyr::transmute(variant, ensg = ensgs, ...)
  }

  # missense variants
  gxv$missense <- variants %>% intersect_REVEL(missense, value = score)

  # nonsense variants
  gxv$nonsense <- variants %>% intersect_REVEL(nonsense)

  # splice site variants
  gxv$splicesite <- variants %>% intersect_REVEL(splicesite)

  # prefix names
  names(gxv) <- paste0("gxv_", names(gxv))
  return(gxv)
}
