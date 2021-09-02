#' Fisher test of significance of overlaps between grouped sets of intervals
#'
#' Calculate Fisher's test on the number of intervals that are shared and unique between two input tibble BEDs.
#' If a grouping variable is provided for bedA, Fisher tests can be performed within these groups of intervals in bedA.
#'
#' This function is built as a modification of the valr::bed_fisher() function.
#' Equivalent to bedtools fisher -a A.bed -b B.bed (performed within each variable group for A.bed)
#'
#' @param bedA A tibble BED with columns: chrom, start, end, ...(metadata cols)...
#' @param bedB A tibble BED with columns: chrom, start, end, ...(metadata cols)...
#' @param bedA_groups A character vector of data variable(s) in bedA upon which to group data and perform the test on each group of intervals separately
#' @param genome A tibble BED of chromosome sizes with columns: chrom, size
#' @param ... Expressions that return a logical value upon which to filter the Fisher results table (thresholds for test output values)
#'
#' @return The outcome of the Fisher test(s) on the overlaps between bedA and bedB
#' @export
bed_fisher_grouped <- function(bedA, bedB, bedA_groups = NULL, genome, ...){
  A <- bedA %>% dplyr::group_by(dplyr::across(dplyr::all_of(bedA_groups)))
  B <- bedB

  A %>%
    dplyr::group_modify(~ valr::bed_fisher(.x, B, genome)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(...)
}
