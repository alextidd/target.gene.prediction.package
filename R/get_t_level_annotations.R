get_t_level_annotations <- function(TSSs,
                                    enriched){

  t <- list()

  # intersect with H3K27ac
  t <- t %>%
    intersect_H3K27ac(TSSs %>% dplyr::select(chrom:end, enst),
                      enriched$H3K27ac,
                      enst)

  # return
  names(t) <- paste0("t_", names(t))
  return(t)
}
