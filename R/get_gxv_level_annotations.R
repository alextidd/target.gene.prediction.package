get_gxv_level_annotations <- function(txv,
                                      variants,
                                      DHSs_master,
                                      specific_DHSs_closest_specific_genes){

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

  # variants in specific DHSs x genes in specific DHSs
  gxv$gxv_specific_DHSs_closest_specific_genes <- variants %>%
    target.gene.prediction.package::bed_intersect_left(DHSs_master, keepBcoords = F) %>%
    dplyr::inner_join(specific_DHSs_closest_specific_genes) %>%
    tidyr::pivot_longer(-c(chrom:DHS), names_to = "celltype", values_to = "symbol") %>%
    dplyr::filter(!is.na(symbol)) %>%
    dplyr::select(variant, celltype, symbol) %>%
    dplyr::distinct() %>%
    dplyr::mutate(value = 1) %>%
    tidyr::pivot_wider(id_cols = -c(celltype, value),
                       names_from = celltype, values_from = value)
  return(gxv)
}
