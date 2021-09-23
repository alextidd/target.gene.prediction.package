get_c_level_annotations <- function(open.variants = open_variants) {
  cat("Annotating credible sets...\n")

  c <- list()

  c$annotations <- open.variants %>%
    dplyr::group_by(cs) %>%
    dplyr::transmute(annotation.name = "inverse_n_variants_within_CS",
                     annotation.value = 1/dplyr::n_distinct(variant),
                     annotation.weight = 1) %>%
    dplyr::distinct()

  # return
  return(c)
}
