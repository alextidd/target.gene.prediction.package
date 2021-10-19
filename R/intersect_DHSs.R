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
intersect_DHSs <- function(l,
                           query,
                           DHSs,
                           ...){

  # intersect with query bed
  l <- DHSs %>%
    lapply(function(x){
      score_cols <- setdiff(colnames(x), c("chrom", "start", "end", "DHS"))
      target.gene.prediction.package::bed_intersect_left(
          bedA = query,
          bedB = x %>% dplyr::select(-DHS),
          keepBcoords = F) %>%
        dplyr::select(-c(chrom:end))%>%
        # group by query column(s) (unit(s) of the annotation)
        dplyr::group_by(...) %>%
        # for overlapping - get max bin to avoid duplicates
        dplyr::summarise(dplyr::across(score_cols, max)) %>%
        dplyr::ungroup()
    })
  names(l) <- paste0("DHSs_", names(l))
  return(l)

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
