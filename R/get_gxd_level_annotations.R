get_gxd_level_annotations <- function(variants,
                                      DHSs_master,
                                      specific_DHSs_closest_specific_genes){
  gxd <- list()
  gxd[["gxd_specific_DHSs_closest_specific_genes"]] <- variants %>%
    bed_intersect_left(DHSs_master, keepBcoords = F) %>%
    dplyr::distinct(DHS) %>%
    dplyr::inner_join(specific_DHSs_closest_specific_genes) %>%
    tidyr::pivot_longer(-DHS, names_to = "celltype", values_to = "symbol") %>%
    dplyr::filter(!is.na(symbol)) %>%
    dplyr::mutate(value = 1) %>%
    tidyr::pivot_wider(id_cols = -c(celltype, value), names_from = celltype, values_from = value)
  return(gxd)
}
