matricise_by_pair <- function(df,
                              txv_master){
  df %>%
    dplyr::ungroup() %>%
    dplyr::right_join(txv_master) %>%
    dplyr::distinct(dplyr::across(c("pair", setdiff(colnames(df), c("DHS",names(txv_master)))))) %>%
    dplyr::mutate_if(is.numeric, tidyr::replace_na, replace = 0) %>%
    tibble::column_to_rownames("pair") %>%
    as.matrix()

}
