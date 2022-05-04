get_vxg_level_annotations <- function(variants,
                                      vxt_master,
                                      enriched){

  vxg <- list()

  # variant-gene distance
  # (between the variant and all of the gene's transcripts' TSSs)
  distance <- vxt_master %>%
    dplyr::group_by(cs, variant, symbol) %>%
    dplyr::summarise(distance = min(distance)) %>%
    dplyr::group_by(variant) %>%
    dplyr::transmute(cs, variant, symbol,
                     inv_distance = dplyr::case_when(distance == 0 ~ 1,
                                                     TRUE ~ 1 / distance),
                     inv_distance_rank = 1 / rank(distance, ties.method = "min")) %>%
    dplyr::ungroup() %>%
    dplyr::distinct()

  vxg$inv_distance <- distance %>%
    dplyr::transmute(cs, variant, symbol,
                     value = inv_distance)

  vxg$inv_distance_rank <- distance %>%
    dplyr::transmute(cs, variant, symbol,
                     value = inv_distance_rank)

  # REVEL #
  intersect_REVEL <- function(variants, revel_df, ...){
    variants %>%
      dplyr::inner_join(revel_df, by = c("chrom" = "chrom", "end" = "position")) %>%
      tidyr::separate_rows(ensgs) %>%
      dplyr::filter(...) %>%
      dplyr::transmute(cs, variant, ensg = ensgs)
  }

  # missense variants
  vxg$missense <- variants %>% intersect_REVEL(missense, score > 0.5)

  # nonsense variants
  vxg$nonsense <- variants %>% intersect_REVEL(nonsense)

  # splice site variants
  vxg$splicesite <- variants %>% intersect_REVEL(splicesite)

  # prefix names
  names(vxg) <- paste0("vxg_", names(vxg))
  return(vxg)
}
