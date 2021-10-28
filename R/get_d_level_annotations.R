get_d_level_annotations <- function(variants){
  cat("Annotating DHSs...\n")

  d <- list()

  d$inv_n_variants <- variants %>%
    bed_intersect_left(DHSs_master, keepBcoords = F) %>%
    dplyr::group_by(DHS) %>%
    dplyr::summarise(value = 1/dplyr::n_distinct(variant))

  # return
  names(d) <- paste0("d_", names(d))
  return(d)
}
