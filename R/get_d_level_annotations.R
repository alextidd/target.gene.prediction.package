get_d_level_annotations <- function(open.variants = open_variants){
  cat("Annotating DHSs...\n")

  d <- list()

  d$annotations <- open.variants %>%
    dplyr::group_by(DHS) %>%
    dplyr::transmute(DHS,
                     annotation.name = "inverse_n_variants_within_DHS",
                     annotation.value = 1/dplyr::n_distinct(variant),
                     annotation.weight = 1) %>%
    dplyr::distinct()

  # return
  return(d)
}
