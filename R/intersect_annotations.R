#' Intersect with the DHS annotations
#'
#' Intersect a query BED with the DHS annotations.
#'
#' @param query A query bed
#' @param annotations A list of BED tibbles to be intersected by the query bed (metadata columns will form the annotation name)
#'
#' @return A tibble of intersected columns, with .query and .annotation suffixes
#' @export
intersect_annotations <- function(query, annots = annotations){

  for(annotation in names(annots)){
    cat("Intersecting query BED with", annotation, "\n")
      df <- annots[[annotation]] %>%
        # add annotation type column
        dplyr::mutate(type = annotation) %>%
        # unite annotations columns together into a single column (for the master table)
        tidyr::unite("annotation", c(type, setdiff(names(.), c("chrom", "start", "end")))) %>%
        # intersect with query bed
        valr::bed_intersect(query, .,
                            suffix = c(".query", ".annotation")) %>%
        dplyr::distinct()
  }

  return(df)
}
