get_cxt_level_annotations <- function(vxt,
                                      variants) {

  # multiHiChIP statistics within each gene-x-cs-x-experiment combination
  multiHiChIP <- vxt[grepl("HiChIP", names(vxt)) & !grepl("binary", names(vxt))] %>%
    purrr::map(~ dplyr::group_by(., cs, enst))

  # count number of loops between gene and CS
  n_multi <- multiHiChIP %>%
    purrr::map(~ dplyr::summarise(., dplyr::across(where(is.numeric), ~ sum(!is.na(.x)))))
  names(n_multi) <- names(n_multi) %>% gsub("vxt_", "n_multi", .)

  # binary - multiHiChIP y/n
  binary_multi <- n_multi %>%
    purrr::map(~ dplyr::mutate(., dplyr::across(where(is.numeric), ~ as.numeric(.x > 1))))
  names(binary_multi) <- names(binary_multi) %>% gsub("multiHiChIP", "multiHiChIP_binary", .)

  # sum of loop values between gene and CS
  sum_multi <- multiHiChIP %>%
    purrr::map(~ dplyr::summarise(., dplyr::across(where(is.numeric), ~ sum(.x, na.rm = T))))
  names(sum_multi) <- names(sum_multi) %>% gsub("vxt_", "sum_multi", .)

  cxt <- c(n_multi,
           binary_multi,
           sum_multi)

  # return
  names(cxt) <- paste0("cxt_", names(cxt))
  return(cxt)
}
