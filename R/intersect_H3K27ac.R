#' Intersect with the DHS annotations
#'
#' Intersect a query BED with the DHS annotations.
#'
#' @param query A query bed
#' @param H3K27ac A BED tibble (the H3K27ac-in-DHSs annotations) to be intersected by the query bed
#' @param ... query columns to be retained in the output
#'
#' @return A tibble of intersected columns, with .query and .annotation suffixes
#' @export
intersect_H3K27ac <- function(l,
                              query,
                              H3K27ac,
                              ...){

  # intersect with query bed
  l <- H3K27ac %>%
    lapply(function(x){
      score_cols <- setdiff(colnames(x), c("chrom", "start", "end", "DHS"))
      bed_intersect_left(
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
  names(l) <- paste0("H3K27ac_", names(l))
  return(l)

}

