get_gxd_level_annotations <- function(variants,
                                      DHSs,
                                      specific_H3K27ac_closest_specific_genes){
  gxd <- list()
  gxd$specific_H3K27ac_closest_specific_genes <- variants %>%
    bed_intersect_left(DHSs, keepBcoords = F) %>%
    dplyr::distinct(DHS) %>%
    dplyr::inner_join(specific_H3K27ac_closest_specific_genes) %>%
    tidyr::pivot_longer(-DHS, names_to = "celltype", values_to = "symbol") %>%
    dplyr::filter(!is.na(symbol)) %>%
    tidyr::pivot_wider(id_cols = -c(celltype, value), names_from = celltype, values_from = value)

  # return
  names(gxd) <- paste0("gxd_", names(gxd))
  return(gxd)
}
