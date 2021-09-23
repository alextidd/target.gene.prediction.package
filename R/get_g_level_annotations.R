get_g_level_annotations <- function(enriched.DHSs = enriched$DHSs){
  cat("Annotating genes...\n")

  g <- list()

  # intersect with DHSs
  g$annotations <- target.gene.prediction.package::intersect_DHSs(
      query = target.gene.prediction.package::TSSs,
      DHSs = enriched.DHSs %>% dplyr::select(chrom:DHS, dplyr::starts_with("signal")),
      enst) %>%
    dplyr::mutate(annotation.weight = 1)

  # return
  return(g)
}
