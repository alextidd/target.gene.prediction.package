get_c_level_annotations <- function(variants) {
  cat("Annotating credible sets...\n")

  c <- list()

  c$inv_n_variants <- variants %>%
    dplyr::group_by(cs) %>%
    dplyr::summarise(value = 1/dplyr::n_distinct(variant))

  # return
  names(c) <- paste0("c_", names(c))
  return(c)
}
