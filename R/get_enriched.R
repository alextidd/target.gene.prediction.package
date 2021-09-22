get_enriched <- function(DHSs, DHSs_metadata, variants){
  enriched <- list()
  enriched[["celltypes"]] <- DHSs %>%
    # Fisher enrichment test of variants in upper-quartile cell-type-specificic H3K27ac marks in DHSs
    dplyr::filter(Method == "specificity",
                  Mark == "H3K27ac",
                  decile == 10) %>%
    target.gene.prediction.package::bed_fisher_grouped(
      bedA = .,
      bedA_groups = "CellType",
      bedB = variants,
      genome = target.gene.prediction.package::ChrSizes,
      # filter for effect and significance
      estimate > 2,
      p.value < 0.05
    ) %>%
    # Extract enriched cell types
    dplyr::pull(CellType) %>%
    {dplyr::filter(DHSs_metadata, code %in% .)}
  enriched[["tissues"]] <- DHSs_metadata %>%
    dplyr::filter(tissue %in% enriched$celltypes$tissue)
  return(enriched)
}
