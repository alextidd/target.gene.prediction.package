#' Intersect with the DHS annotations
#'
#' Intersect a query BED with the DHS annotations.
#'
#' @param query A query bed
#' @param DHSs A BED tibble (the DHSs) to be intersected by the query bed
#' @param ... query columns to be retained in the output
#'
#' @return A tibble of intersected columns, with .query and .annotation suffixes
#' @export
intersect_DHSs <- function(query,
                           DHSs,
                           ...){

  cat("Intersecting query BED with DHSs\n")
  # intersect with query bed
  target.gene.prediction.package::bed_intersect_left(
    query,
    DHSs,
    keepBcoords = F) %>%
    tidyr::gather(key = "annotation.name",
                  value = "annotation.value",
                  -c(names(query), "DHS")) %>%
    dplyr::transmute(...,
                     annotation.name = paste0("DHSs_", annotation.name),
                     annotation.value)
}
