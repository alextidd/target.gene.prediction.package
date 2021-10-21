get_txc_level_annotations <- function(txv,
                                      variants) {
  cat("Annotating transcript x credible set pairs...\n")

  # multicontact statistics within each gene-x-cs-x-experiment combination
  multicontact <- txv[grepl("contact", names(txv)) & !grepl("binary", names(txv))] %>%
    purrr::map(~ dplyr::left_join(., variants %>% dplyr::select(variant, cs)) %>%
                 dplyr::group_by(cs, enst))

  # count number of loops between gene and CS
  n_multi <- multicontact %>%
    purrr::map(~ dplyr::summarise(., dplyr::across(where(is.numeric), ~ sum(!is.na(.x)))))
  names(n_multi) <- names(n_multi) %>% gsub("txv_", "txc_n_multi", .)

  # binary - multicontact y/n
  binary_multi <- n_multi %>%
    purrr::map(~ dplyr::mutate(., dplyr::across(where(is.numeric), ~ as.numeric(.x > 1))))
  names(binary_multi) <- names(binary_multi) %>% gsub("multicontact", "multicontact_binary", .)

  # sum of loop values between gene and CS
  sum_multi <- multicontact %>%
    purrr::map(~ dplyr::summarise(., dplyr::across(where(is.numeric), ~ sum(.x, na.rm = T))))
  names(sum_multi) <- names(sum_multi) %>% gsub("txv_", "txc_sum_multi", .)

  txc <- c(n_multi,
           binary_multi,
           sum_multi)

  # return
  return(txc)
}
