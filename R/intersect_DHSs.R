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
intersect_DHSs <- function(list,
                           query,
                           DHSs,
                           ...){

  # intersect with query bed
  list[["DHSs_specificity"]] <- int_func(query, DHSs[["specificity"]])
  list[["DHSs_signal"]] <- int_func(query, DHSs[["signal"]])
  return(list)

}

# internal generic intersection function
int_func <- function(query, DHSs){
  target.gene.prediction.package::bed_intersect_left(
    query,
    DHSs %>% dplyr::select(-DHS),
    keepBcoords = F) %>%
    dplyr::select(-c(chrom, start, end))
}


# target.gene.prediction.package::bed_intersect_left(
#   query,
#   DHSs,
#   keepBcoords = F) %>%
#   tidyr::gather(key = "annotation.name",
#                 value = "annotation.value",
#                 -c(names(query), "DHS")) %>%
#   dplyr::transmute(...,
#                    annotation.name = paste0("DHSs_", annotation.name),
#                    annotation.value)
