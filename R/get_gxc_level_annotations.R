get_gxc_level_annotations <- function() {
  cat("Annotating gene x", trait, "credible set pairs...\n")

  gxc <- list()

  multicontact <- gxv$contact %>%
    # multicontact statistics within each gene-x-cs-x-experiment combination
    dplyr::group_by(cs, enst, annotation.name) %>%
    dplyr::mutate(
      # number of gxc loops within each experiment
      inv_n_contacts = 1/dplyr::n_distinct(InteractionID),
      # sum of the scores of gxc loops within each experiment
      inv_sum_contacts = 1/sum(annotation.value),
      # collapse Interaction IDs
      InteractionID = paste0(unique(InteractionID), collapse = ",")
    )

  gxc$n_multicontact <- multicontact %>%
    dplyr::ungroup()  %>%
    dplyr::transmute(cs,
                     enst,
                     InteractionID,
                     annotation.name = paste0("n_multicontact_", annotation.name),
                     annotation.value = inv_n_contacts,
                     # weighting (more for enriched tissues)
                     annotation.weight = dplyr::case_when(annotation.name %in% enriched_contact_elements ~ 2 * weight$gxc$multicontact,
                                                          TRUE ~ weight$gxc$multicontact)) %>%
    dplyr::distinct()

  gxc$sum_multicontact <- multicontact %>%
    dplyr::ungroup()  %>%
    dplyr::transmute(cs,
                     enst,
                     InteractionID,
                     annotation.name = paste0("sum_multicontact_", annotation.name),
                     annotation.value = inv_sum_contacts,
                     # weighting (more for enriched tissues)
                     annotation.weight = dplyr::case_when(annotation.name %in% enriched_contact_elements ~ 2 * weight$gxc$multicontact,
                                                          TRUE ~ weight$gxc$multicontact)) %>%
    dplyr::distinct()

  # return
  return(gxc)
}
