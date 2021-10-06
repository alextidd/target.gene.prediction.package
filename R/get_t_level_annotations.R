get_t_level_annotations <- function(enriched){
  cat("Annotating transcripts...\n")

  t <- list()

  # intersect with DHSs
  t <- t %>%
    intersect_DHSs(target.gene.prediction.package::TSSs %>% dplyr::select(chrom:end, enst),
                   enriched$DHSs)

  # return
  names(t) <- paste0("t_", names(t))
  return(t)
}
