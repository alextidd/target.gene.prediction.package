#' Left join intersection of two tibble BEDs
#'
#' Intersect two tibble BEDs and retain the full intervals from the lefthand tibble BED that have any overlap.
#' If the two input tibble BEDs share any common metadata column names, suffixes must be provided.
#' This function is built as a modification of the valr::bed_intersect() function.
#' Equivalent to bedtools -a A -b B -wa.
#'
#' @param bedA A tibble BED with columns: chrom, start, end, ...(metadata cols)...
#' @param bedB A tibble BED with columns: chrom, start, end, ...(metadata cols)...
#' @param suffix A character vector of suffixes to add to metadata column names. (If the two input tibble BEDs share any common metadata column names, suffixes must be provided.)
#'
#' @return A tibble BED of all intervals from bedA that contain intersects with bedB, plus all metadata columns from both
#' @export
bed_intersect_left <- function(bedA, bedB, suffix = c("","")){

  # Make sure suffix has been provided if there are common metadata columns
  check_A_vs_B_metadata_cols(bedA, bedB)

  # Left join intersect
  bedA %>%
    valr::bed_intersect(bedB) %>%
    dplyr::select(-c("start.y", "end.y")) %>%
    dplyr::rename_with(~ gsub("[.]x", suffix[1], .)) %>%
    dplyr::rename_with(~ gsub("[.]y", suffix[2], .))

}

check_A_vs_B_metadata_cols <- function(bedA, bedB){
  suffix <- NULL
  metaA <- names(bedA) %>% setdiff(c("chrom","start","end"))
  metaB <- names(bedB) %>% setdiff(c("chrom","start","end"))
  if(length(intersect(metaA, metaB) > 0) && (suffix == c("",""))){
    stop("bedA and bedB have common metadata column name(s) but suffixes have not been provided.",
         "\nPlease provide suffixes to be appended to the metadata columns in the output.",
         "\nShared metadata column name(s): ", intersect(metaA, metaB))
  }
  if(length(intersect(metaA, metaB) > 0) && suffix[1] == suffix[2]){
    stop("bedA and bedB have common metadata column name(s) but unique suffixes have not been provided.",
         "\nPlease provide two unique suffixes to be appended to the metadata columns in the output.",
         "\nSuffixes: ", suffix)
  }
}
