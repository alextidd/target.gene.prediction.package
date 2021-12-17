matricise_by_pair <- function(df,
                              txv_master){
  # add value = 1 if there is no value (numeric) column
  if(df %>% dplyr::ungroup() %>% dplyr::select(where(is.numeric)) %>% ncol == 0){df$value <- 1}

  # MA needs all rows in the same order
  mat <- df %>%
    dplyr::rowwise() %>%
    # aggregate per pair - mean across samples
    dplyr::mutate(value = mean(dplyr::c_across(where(is.numeric)), na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::right_join(txv_master, by = intersect(names(txv_master), names(df))) %>%
    dplyr::select(pair, value) %>%
    dplyr::mutate_if(is.numeric, tidyr::replace_na, replace = 0) %>%
    tibble::column_to_rownames("pair") %>%
    as.matrix

  return(mat)
}
