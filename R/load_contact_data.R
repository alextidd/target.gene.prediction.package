#' Load contact (3C) data
#'
#' Load the directory of contact data, downloaded from an FTP(?) alongside the package.
#' The function will take all files in this directory ending in *.bedpe, convert them to
#' BEDPE list objects, add an InteractionID column and normalise the score column within
#' each sample.
#'
#' @param contactDir The directory in which the package's accompanying contact data is stored.
#' Files must have columns ( chrA | startA | endA | chrB | startB | endB | score | celltype ) and *.bedpe file extension.
#'
#' @return `contact` object (list of BEDPE lists containing loops and their normalised scores for multiple cell types)
#' @export
load_contact_data <- function(contactDir){
  contact <- list()
  for(file in list.files(contactDir, pattern = "bedpe", full.names = T)){
    info <- basename(file) %>%  sub("(.*)\\..*$", "\\1", .)
    contact[[info]] <- target.gene.prediction.package::import_BEDPE_to_List(
      file, metadata_cols = c("score", "CellType")
      ) %>%
      purrr::map(~ .x %>%
                   # for infinite score values, set equal to the maximum non-infinite score
                   dplyr::mutate(score = dplyr::case_when(is.infinite(score) ~ max(score[!is.infinite(score)]),
                                                          TRUE ~ score)) %>%
                   # normalise score col to maximum
                   dplyr::mutate(score = score/max(score)))
  }
  return(contact)
}
