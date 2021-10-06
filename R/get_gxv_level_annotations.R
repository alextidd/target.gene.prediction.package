get_gxv_level_annotations <- function(txv){

  cat("Annotating gene x variant pairs...\n")

  gxv <- list()

  gxv$gxv_inv_distance_rank <- txv$txv_inv_distance %>%
    dplyr::left_join(target.gene.prediction.package::TSSs %>% dplyr::select(enst, symbol)) %>%
    # Rank for each variant-gene pair
    dplyr::group_by(symbol) %>%
    dplyr::mutate(closest_distance = max(value)) %>%
    dplyr::distinct(variant, symbol, closest_distance) %>%
    dplyr::group_by(variant) %>%
    dplyr::transmute(
      symbol,
      variant,
      value = 1/rank(closest_distance, ties.method = "min")
      )

  return(gxv)
}
