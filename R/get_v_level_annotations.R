get_v_level_annotations <- function(variants,
                                    enriched,
                                    txv_master){

  v <- list()

  # intersect with DHS binnings
  v <- v %>%
    intersect_DHSs(query = variants,
                   DHSs = enriched$DHSs,
                   variant)

  # calculate n genes near each variant
  v <- txv_master %>%
    dplyr::group_by(variant) %>%
    dplyr::summarise(inv_n_genes = 1/dplyr::n_distinct(symbol),
                     inv_n_transcripts = 1/dplyr::n_distinct(enst)) %>%
    tidyr::pivot_longer(c(inv_n_genes, inv_n_transcripts),
                        names_to = "split",
                        values_to = "value") %>%
    split(as.factor(.$split)) %>%
    purrr::map(~ dplyr::select(., -split)) %>%
    c(., v)

  # return
  names(v) <- paste0("v_", names(v))
  return(v)
}
