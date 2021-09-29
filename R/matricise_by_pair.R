matricise_by_pair <- function(df,
                              gxv_master){
  df %>%
    dplyr::ungroup() %>%
    dplyr::right_join(gxv_master) %>%
    dplyr::select(c("pair", setdiff(colnames(df), names(gxv_master)))) %>%
    dplyr::mutate_if(is.numeric, tidyr::replace_na, replace = 0) %>%
    tibble::column_to_rownames("pair") %>%
    as.matrix()
}
