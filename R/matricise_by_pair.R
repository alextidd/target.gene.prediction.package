matricise_by_pair <- function(df,
                              txv_master){
  # MA needs all rows in the same order
  dplyr::left_join(txv_master, dplyr::ungroup(df), by = intersect(names(txv_master), names(df))) %>%
    dplyr::select(c("pair", setdiff(colnames(df), c("DHS",names(txv_master))))) %>%
    dplyr::mutate_if(is.numeric, tidyr::replace_na, replace = 0) %>%
    tibble::column_to_rownames("pair") %>%
    as.matrix()

}
