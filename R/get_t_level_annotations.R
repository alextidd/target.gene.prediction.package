get_t_level_annotations <- function(TSSs,
                                    enriched_DHSs){
  cat("Annotating transcripts...\n")

  t <- list()

  # intersect with DHSs
  t <- t %>%
    intersect_DHSs(TSSs %>% dplyr::select(chrom:end, enst),
                   enriched_DHSs,
                   enst)

  # return
  names(t) <- paste0("t_", names(t))
  return(t)
}
