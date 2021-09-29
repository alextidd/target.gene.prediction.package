get_d_level_annotations <- function(open.variants = open_variants){
  cat("Annotating DHSs...\n")

  d <- list()

  d$inv_n_variants <- open.variants %>%
    dplyr::group_by(DHS) %>%
    dplyr::transmute(DHS,
                     value = 1/dplyr::n_distinct(variant)) %>%
    dplyr::distinct()


  # return
  names(d) <- paste0("d_", names(d))
  return(d)
}
