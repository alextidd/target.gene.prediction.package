#' Bind and widen annotation tibbles
#'
#' Bind annotation tibbles and then apply pivot_wider to the annotation columns, converting from long-format
#' to wide-format, with only 1 variant or gene or variant-gene pair per row and 1 annotation per column
#'
#' @param id_cols Vector of column names that uniquely identifies each observation, sent to pivot_wider(id_cols)
#' @param annotation.level Name of the annotation level (e.g. pair, gene, variant) that the annotations fall within
#' @param ... The long-format annotation tibble(s) to bind and widen
#'
#' @return Wide-format annotation dataframe, with only one observation per row and only one annotation per column
#' @export
bind_and_widen_annotations <- function(id_cols, annotation_level, ...){

  dplyr::bind_rows(...) %>%
    dplyr::mutate(annotation.level = annotation_level) %>%
    tidyr::pivot_wider(id_cols = id_cols,
                       names_from = c(annotation.level, annotation.name),
                       values_from = annotation.value) %>%
    readr::type_convert() ## TODO: silence column specification output

}
