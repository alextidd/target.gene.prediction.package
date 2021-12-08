matricise_by_pair <- function(df,
                              txv_master){
  # MA needs all rows in the same order
  df %>%
    dplyr::rowwise() %>%
    # aggregate per pair - maximum across samples
    dplyr::mutate(value = max(dplyr::c_across(where(is.numeric)), na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::right_join(txv_master, by = intersect(names(txv_master), names(df))) %>%
    dplyr::select(pair, value) %>%
    dplyr::mutate_if(is.numeric, tidyr::replace_na, replace = 0) %>%
    tibble::column_to_rownames("pair") %>%
    as.matrix
}
