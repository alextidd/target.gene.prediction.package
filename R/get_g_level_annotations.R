get_g_level_annotations <- function(txv_master,
                                    enriched){
  cat("Annotating transcripts...\n")

  g <- list()

  # intersect with DHSs
  g$g_expression_binary <- enriched$expression %>%
    tibble::as_tibble(rownames = "ensg") %>%
    dplyr::inner_join(txv_master %>% dplyr::distinct(ensg))

  # return
  return(g)
}
