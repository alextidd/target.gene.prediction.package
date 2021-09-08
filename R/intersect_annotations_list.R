#' Intersect with the annotations list
#'
#' Intersect a query BED with the annotations list.
#'
#' @param query A query bed
#'
#' @return A tibble of intersected columns, with .query and .annotation suffixes
#' @export
intersect_annotations_list <- function(query){
  df <- list()
  for(annotation in names(target.gene.prediction.package::annotations)){
    cat("============================================\n",
        "Intersecting query BED with", annotation, "\n")
    df[[annotation]] <- valr::bed_intersect(query,
                                            target.gene.prediction.package::annotations[[annotation]],
                                            suffix = c(".query", ".annotation"))
  }
  dplyr::bind_rows(df)
}
