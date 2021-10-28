get_txd_level_annotations <- function(ensts_near_vars,
                                      DHSs_master) {
  cat("Annotating transcript x DHS pairs...\n")

  txd <- list()

  #Wrong - TODO must fix
  # closest DHS to the transcript
  txd$closest <- valr::bed_closest(TSSs %>%
                                     dplyr::filter(enst %in% ensts_near_vars),
                                   DHSs_master) %>%
    # valr's .overlap and .dist results are unreliable, generating my own
    dplyr::rowwise() %>%
    dplyr::mutate(dist = min(abs(end.y - start.x), abs(end.x - start.y)),
                  overlap = start.x > start.y & end.x < end.y) %>%
    # get closest / overlapping
    dplyr::group_by(enst.x) %>%
    dplyr::filter(overlap | overlap & dist == min(dist) | dist == min(dist)) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(DHS = DHS.y,
                     enst = enst.x,
                     value = 1)

  # return
  names(txd) <- paste0("txd_", names(txd))
  return(txd)
}
