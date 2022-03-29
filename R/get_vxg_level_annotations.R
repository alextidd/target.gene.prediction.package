get_vxg_level_annotations <- function(variants,
                                      vxt_master,
                                      enriched){

  vxg <- list()

  # variant-gene closest distance (among all of the gene's transcripts' TSSs)
  vxg$inv_distance_rank <- vxt_master %>%
    dplyr::group_by(cs, variant, symbol) %>%
    dplyr::summarise(distance = min(distance)) %>%
    dplyr::group_by(variant) %>%
    dplyr::transmute(
      cs, variant, symbol,
      value = 1/rank(distance, ties.method = "min")
    )

  # REVEL #
  intersect_REVEL <- function(variants, revel_df, ...){
    variants %>%
      dplyr::inner_join(revel_df, by = c("chrom" = "chrom", "end" = "position")) %>%
      tidyr::separate_rows(ensgs) %>%
      dplyr::transmute(cs, variant, ensg = ensgs, ...)
  }

  # missense variants
  vxg$missense <- variants %>% intersect_REVEL(missense, value = score)

  # nonsense variants
  vxg$nonsense <- variants %>% intersect_REVEL(nonsense)

  # splice site variants
  vxg$splicesite <- variants %>% intersect_REVEL(splicesite)

  # prefix names
  names(vxg) <- paste0("vxg_", names(vxg))
  return(vxg)
}
