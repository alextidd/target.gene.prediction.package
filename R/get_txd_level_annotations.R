get_txd_level_annotations <- function(ensts_near_vars,
                                      DHSs_master) {
  cat("Annotating transcript x DHS pairs...\n")

  txd <- list()

  # closest DHS to the transcript
  txd$closest <- valr::bed_closest(target.gene.prediction.package::TSSs %>%
                                     dplyr::filter(enst %in% ensts_near_vars),
                                   DHSs_master) %>%
    dplyr::group_by(enst.x) %>%
    dplyr::filter(abs(.dist) == min(abs(.dist))) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(DHS = DHS.y,
                     enst = enst.x,
                     value = 1)

  # return
  names(txd) <- paste0("txd_", names(txd))
  return(txd)
}
