get_gxd_level_annotations <- function(gxv.contact = gxv$contact,
                                      ensts.near.vars = gxv$distance_score$enst,
                                      weight.gxd = weight$gxd,
                                      .DHSs = DHSs) {
  cat("Annotating gene x DHS pairs...\n")

  gxd <- list()

  multicontact <- gxv.contact %>%
    # multicontact statistics within each gene-x-dhs-x-experiment combination
    dplyr::group_by(DHS, enst, annotation.name) %>%
    dplyr::mutate(
      #  count of gxd loops
      inv_n_contacts = 1/dplyr::n_distinct(InteractionID),
      # sum of gxd loop values
      inv_sum_contacts = 1/sum(annotation.value),
      # collapse Interaction IDs
      InteractionID = paste0(unique(InteractionID), collapse = ",")
    )

  gxd$n_multicontact <- multicontact %>%
    dplyr::ungroup()  %>%
    dplyr::transmute(DHS,
                     enst,
                     InteractionID,
                     annotation.name = paste0("n_multicontact_", annotation.name),
                     annotation.value = inv_n_contacts,
                     annotation.weight = weight.gxd$multicontact) %>%
    dplyr::distinct()

  gxd$sum_multicontact <- multicontact %>%
    dplyr::ungroup()  %>%
    dplyr::transmute(DHS,
                     enst,
                     InteractionID,
                     annotation.name = paste0("sum_multicontact_", annotation.name),
                     annotation.value = inv_sum_contacts,
                     # weighting (more for enriched tissues)
                     annotation.weight = weight.gxd$multicontact) %>%
    dplyr::distinct()

  # Closest DHS to the gene?
  gxd$closest <- valr::bed_closest(target.gene.prediction.package::TSSs %>%
                                     dplyr::filter(enst %in% ensts.near.vars),
                                   .DHSs %>%
                                     dplyr::select(chrom:end, DHS) %>%
                                     dplyr::distinct()
                                   ) %>%
    dplyr::group_by(enst.x) %>%
    dplyr::filter(abs(.dist) == min(abs(.dist))) %>%
    dplyr::select(DHS = DHS.y,
                  enst = enst.x) %>%
    dplyr::distinct() %>%
    dplyr::transmute(DHS, enst,
                     annotation.name = "closest",
                     annotation.value = 1,
                     annotation.weight = weight.gxd$closest)

  # return
  return(gxd)
}
