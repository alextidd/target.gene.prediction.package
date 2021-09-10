#' Bind and weight and widen annotation tibbles
#'
#' Bind annotation tibbles, optionally weight values by a weighting column, and then apply pivot_wider
#' to the annotation columns, converting from long-format to wide-format, with only 1 variant or gene or
#' variant-gene pair per row and 1 annotation per column
#'
#' @param id_cols Vector of column names that uniquely identifies each observation, sent to pivot_wider(id_cols)
#' @param annotation_level Name of the annotation level (e.g. pair, gene, variant) that the annotations fall within
#' @param ... The long-format annotation tibble(s) to bind and widen
#'
#' @return Wide-format annotation dataframe, with only one observation per row and only one annotation per column
#' @export
bind_and_weight_and_widen_annotations <- function(id_cols, annotation_level, ...){

  # BIND
  dplyr::bind_rows(...) %>%
  # WEIGHT
    dplyr::mutate(annotation.level = annotation_level,
                  annotation.value = annotation.value * annotation.weight) %>%
  # WIDEN
    tidyr::pivot_wider(df,
                       id_cols = id_cols,
                       names_from = c(annotation.level, annotation.name),
                       values_from = annotation.value) %>%
    readr::type_convert() ## TODO: silence column specification output

}




### weighting option - did not work for some reason (default apply_weight=T did not work)
# # Make sure weighting provided if apply_weight = T
# if(apply_weight == TRUE && "annotation.weight" %ni% names(df)) {
#   stop("apply_weight option set to TRUE but no `annotation.weight` column exists in the given annotation tibble(s).\n",
#        "Set apply_weight=FALSE or add/name an `annotation.weight` column in the input (...).")
# }
#
# # WEIGHT
# if(apply_weight == TRUE) {
#   df <- df %>%
#     dplyr::mutate(annotation.value = annotation.value * annotation.weight)
# }

