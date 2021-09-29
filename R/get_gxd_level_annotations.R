get_gxd_level_annotations <- function(ensts.near.vars = unique(gxv$inv_distance$enst),
                                      DHSs.master = DHSs_master) {
  cat("Annotating gene x DHS pairs...\n")

  gxd <- list()

  # Closest DHS to the gene?
  gxd$closest <- valr::bed_closest(target.gene.prediction.package::TSSs %>%
                                     dplyr::filter(enst %in% ensts.near.vars),
                                   DHSs.master) %>%
    dplyr::group_by(enst.x) %>%
    dplyr::filter(abs(.dist) == min(abs(.dist))) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(DHS = DHS.y,
                     enst = enst.x,
                     value = 1)

  # return
  names(gxd) <- paste0("gxd_", names(gxd))
  return(gxd)
}
