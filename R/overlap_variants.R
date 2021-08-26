#' Overlap variants with another BED file
#'
#' Function to overlap the user-provided trait variants file with other relevant BED files from package data
#'
#' @param vars A GRanges object of trait variants with metadata columns snpid and locus
#'
#' @param bed_to_overlap A GRanges object of genomic annotations
#'
#' @return A GRanges object representation of the variant BED file
#' @export
overlap_variants <- function(vars, bed_to_overlap) {
  plyranges::join_overlap_inner(vars, bed_to_overlap)
}
