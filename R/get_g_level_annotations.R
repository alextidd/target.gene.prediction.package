get_g_level_annotations <- function(txv_master,
                                    enriched){

  g <- list()

  # gene signal/spec bins
  g <- enriched$expression %>%
    purrr::map( ~ .x %>% dplyr::inner_join(txv_master %>% dplyr::distinct(ensg), by = "ensg"))
  names(g) <- paste0("expression_", names(g))

  # gene expressed
  g$expressed <- enriched$expressed %>%
    dplyr::inner_join(txv_master %>% dplyr::distinct(ensg), by = "ensg")

  # return
  names(g) <- paste0("g_", names(g))
  return(g)
}
