get_g_level_annotations <- function(enriched.DHSs = enriched$DHSs){
  cat("Annotating genes...\n")

  g <- list()

  # intersect with DHSs
  g <- g %>%
    intersect_DHSs(target.gene.prediction.package::TSSs %>% dplyr::select(chrom:end, enst),
                   enriched.DHSs)

  # return
  names(g) <- paste0("g_", names(g))
  return(g)
}
