get_gxc_level_annotations <- function(.gxv = gxv,
                                      open.variants = open_variants) {
  cat("Annotating gene x credible set pairs...\n")

  # multicontact statistics within each gene-x-cs-x-experiment combination
  multicontact <- .gxv[greplany("contact", names(.gxv))] %>%
    purrr::map(~ dplyr::left_join(., open.variants %>% dplyr::select(variant, cs)) %>%
                 dplyr::group_by(cs, enst))

  # count number of loops between gene and CS
  n_multi <- multicontact %>%
    purrr::map(~ dplyr::summarise(., dplyr::across(where(is.numeric), ~ sum(!is.na(.x)))))
  names(n_multi) <- paste0("n_multi", names(n_multi))

  # sum of loop values between gene and CS
  sum_multi <- multicontact %>%
    purrr::map(~ dplyr::summarise(., dplyr::across(where(is.numeric), ~ sum(.x, na.rm = T))))
  names(sum_multi) <- paste0("sum_multi", names(sum_multi))

  gxc <- c(n_multi,
           sum_multi)

  # return
  names(gxc) <- paste0("gxc_", names(gxc))
  return(gxc)
}
