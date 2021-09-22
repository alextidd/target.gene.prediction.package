get_g_level_annotations <- function(){
  cat("Annotating genes...\n")

  g <- list()

  # intersect with DHSs
  g$annotations <- target.gene.prediction.package::TSSs %>%
    target.gene.prediction.package::intersect_DHSs(
      enst,
      annots = enriched_DHSs %>% dplyr::filter(Method == "signal")) %>%
    dplyr::mutate(annotation.weight = 1)

  # return
  return(g)
}
