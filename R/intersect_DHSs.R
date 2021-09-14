#' Intersect with the DHS annotations
#'
#' Intersect a query BED with the DHS annotations.
#'
#' @param query A query bed
#' @param ... query columns to be retained in the output
#' @param DHSss A BED tibble (the DHSs) to be intersected by the query bed
#' @param annotation.type A name for the type of annotation (`annotation.type` + the metadata columns will form the annotation name)
#'
#' @return A tibble of intersected columns, with .query and .annotation suffixes
#' @export
intersect_DHSs <- function(query,
                           ...,
                           annots = DHSs,
                           specificity.enriched.celltypes = specificity_enriched_celltypes,
                           DHSs.enriched = DHSs_enriched){

  cat("Intersecting query BED with DHSs\n")
  annots %>%
    # add annotation type column
    dplyr::mutate(type = "DHSs") %>%
    # unite annotations columns together into a single column (for the master table)
    tidyr::unite("annotation", c(type, setdiff(names(.), c("chrom", "start", "end"))), remove = F) %>%
    # intersect with query bed
    valr::bed_intersect(query, .,
                        suffix = c("", ".DHSs")) %>%
    dplyr::distinct() %>%
    dplyr::transmute(...,
                     annotation.name = annotation.DHSs,
                     # The value of the annotation is its decile (signal / specificity)
                     annotation.value = decile.DHSs,
                     # Enriched celltype annotations are weighted higher
                     annotation.weight = case_when(
                       # 3 = enriched celltype
                       CellType.DHSs == specificity.enriched.celltypes$code ~ 3,
                       # 2 = celltype from same tissue as enriched celltype
                       CellType.DHSs %in% DHSs.enriched$code ~ 2,
                       # 1 = not enriched
                       TRUE ~ 1))
}
