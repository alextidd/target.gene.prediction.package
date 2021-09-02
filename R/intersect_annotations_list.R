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

    df[[annotation]] <- target.gene.prediction.package::annotations[[annotation]] %>%
      # add annotation type column
      dplyr::mutate(type = annotation) %>%
      # unite annotations columns together into a single column (for the master table)
      tidyr::unite("annotation", c(type, 4:ncol(.))) %>%
      # intersect with query bed
      valr::bed_intersect(query,
                          .,
                          suffix = c(".query", ".annotation"))
  }
  dplyr::bind_rows(df)
}
