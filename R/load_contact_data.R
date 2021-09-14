#' Load contact (3C) data
#'
#' Load the directory of contact data, downloaded from an FTP(?) alongside the package.
#' The function will take all files in this directory ending in *.bedpe, convert them to
#' BEDPE list objects, add an InteractionID column and normalise the score column within
#' each sample.
#'
#' @param referenceDir The directory in which the package's accompanying data is stored
#' (contact data is in the contact/ subdirectory).
#' Files must have columns ( chrA | startA | endA | chrB | startB | endB | score | celltype ) and *.bedpe file extension.
#'
#' @return `contact` object (list of BEDPE lists containing loops and their normalised scores for multiple cell types)
#' @export
load_contact_data <- function(referenceDir){
  contactDir <- paste0(referenceDir, "contact/")
  contact <- list()
  for(file in list.files(contactDir, pattern = "bedpe", full.names = T)){
    info <- basename(file) %>%  sub("(.*)\\..*$", "\\1", .)
    contact[[info]] <- target.gene.prediction.package::import_BEDPE_to_List(
      file, metadata_cols = c("score", "CellType"), prefix_InteractionID = info
      ) %>%
      purrr::map(~ .x %>%
                   # for infinite score values, set equal to the maximum non-infinite score
                   dplyr::mutate(score = dplyr::case_when(is.infinite(score) ~ max(score[!is.infinite(score)]),
                                                          TRUE ~ score)) %>%
                   # normalise score col to maximum
                   dplyr::mutate(score = score/max(score)))

    # find duplicate, mirrored loop entries within each dataset and filter to exclude the non-maximum scores
    IDs_to_exclude <- dplyr::full_join(contact[[info]]$first, contact[[info]]$last, by = "InteractionID") %>%
      dplyr::rowwise()%>%
      dplyr::mutate(loop_ends = paste(min(start.x, start.y), min(end.x, end.y), max(start.x, start.y), max(end.x, end.y), sep = ".")) %>%
      dplyr::group_by(loop_ends) %>%
      dplyr::filter(dplyr::n() > 1,
                    score.x != max(score.x)) %>%
      dplyr::pull(InteractionID)

    if(length(IDs_to_exclude) > 0){
      message(length(IDs_to_exclude), " / ", dplyr::n_distinct(contact[[info]]$first$InteractionID),
              " loops are duplicates. Each group of loops with identical ends is filtered to only include the maximum-scoring loop.")
      contact[[info]] <- contact[[info]] %>%
        purrr::map(~ .x %>%
                     # exclude lower-scoring duplicates
                     dplyr::filter(InteractionID %ni% IDs_to_exclude)
                   )
    }
  }

  return(contact)
}


