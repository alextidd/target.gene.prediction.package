get_g_level_annotations <- function(txv_master,
                                    enriched){
  cat("Annotating transcripts...\n")

  g <- list()

  # gene expressed
  g$expression <- enriched$expression %>%
    tibble::as_tibble(rownames = "ensg") %>%
    dplyr::inner_join(txv_master %>% dplyr::distinct(ensg))

  # return
  names(g) <- paste0("g_", names(g))
  return(g)
}
