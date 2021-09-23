get_gxv_level_annotations <- function(open.variants = open_variants,
                                      variant.to.gene.max.distance = variant_to_gene_max_distance,
                                      weight.gxv = weight$gxv,
                                      enriched.contact_elements = enriched$contact_elements,
                                      .contact = contact) {
  cat("Annotating gene x variant pairs...\n")

  gxv <- list()

  # TSS distance scoring method (these pairings are the only ones to consider)
  distance <- open.variants %>%
    # Get genes within variant.to.gene.max.distance of each variant
    valr::bed_slop(both = variant.to.gene.max.distance,
                   genome = target.gene.prediction.package::ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs,
                        suffix = c(".variant", ".TSS")) %>%
    dplyr::group_by(variant.variant) %>%
    # Calculate inverse of the absolute bp distance for each variant-gene pair
    dplyr::mutate(pair.distance = abs((start.variant + variant.to.gene.max.distance) - start.TSS),
                  pair.inverse_distance = 1/pair.distance,
                  # ranking transcript TSSs (if two transcript TSSs are equidistant to the variant, they will receive the same, lower rank)
                  pair.inverse_distance_rank = 1/rank(pair.distance, ties.method = "min")) %>%
    dplyr::select(variant = variant.variant,
                  cs = cs.variant,
                  DHS = DHS.variant,
                  enst = enst.TSS,
                  pair.inverse_distance,
                  pair.inverse_distance_rank)

  # variant-TSS inverse distance score method
  gxv$distance_score <- distance %>%
    dplyr::transmute(variant, cs, enst, DHS,
                     annotation.name = "inverse_distance_within_max_distance",
                     annotation.value = pair.inverse_distance,
                     annotation.weight = weight.gxv$distance_rank)

  # variant-TSS distance rank method
  gxv$distance_rank <- distance %>%
    dplyr::transmute(variant, cs, enst, DHS,
                     annotation.name = "inverse_distance_rank_within_max_distance",
                     annotation.value = pair.inverse_distance_rank,
                     annotation.weight = weight.gxv$distance_rank)

  # closest variant-TSS method
  gxv$closest <- distance %>%
    dplyr::filter(pair.inverse_distance_rank == 1) %>%
    dplyr::transmute(variant, cs, enst, DHS,
                     annotation.name = "closest",
                     annotation.value = 1,
                     annotation.weight = weight.gxv$closest)

  # intersect loop ends, by cell type, with enhancer variants and gene TSSs
  # (finds interaction loops with a variant at one end and a TSS at the other)
  gxv$contact <- .contact %>%
    purrr::map(~ intersect_BEDPE(
      # ! For mutually exclusive intersection with HiC ranges, make variant intervals 1bp long, equal to the end position
      SNPend = open.variants %>% dplyr::mutate(start = end),
      TSSend = target.gene.prediction.package::TSSs,
      bedpe = .) %>%
        dplyr::transmute(variant, cs, DHS, enst,
                         InteractionID,
                         annotation.value = score)) %>%
    dplyr::bind_rows(.id = "annotation.name") %>%
    dplyr::mutate(
      annotation.weight = dplyr::case_when(annotation.name %in% enriched.contact_elements ~ 2 * weight.gxv$contact,
                                           TRUE ~ weight.gxv$contact)) %>%
    # Make sure all interactions are within 2Mb - hard filter, ignore everything further
    dplyr::inner_join(distance %>% dplyr::select(variant, cs, enst))
  ## ~10% of variant-TSS interactions indicated by the contact data
  ## are further than 2Mb apart and are thus eliminated

  # get variants at promoters - score by sum of signal and specificity per enriched celltype at the DHS (~promoter activity)
  gxv$promoter <- open.variants %>%
    # Get variants within promoter regions
    target.gene.prediction.package::bed_intersect_left(
      target.gene.prediction.package::promoters,
      keepBcoords = F, keepBmetadata = T) %>%
    dplyr::left_join(enriched.DHSs %>% dplyr::select(-c(chrom:end))) %>%
    tidyr::gather(key = "annotation",
                  value = "annotation.value",
                  dplyr::starts_with(c("specificity", "signal"))) %>%
    dplyr::mutate(CellType = annotation %>% gsub(".*_", "", .)) %>%
    dplyr::group_by(variant, enst, CellType) %>%
    dplyr::mutate(annotation.value = sum(annotation.value)) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(variant,
                     enst,
                     annotation.name = paste0(CellType, "_promoter_signal_plus_specificity"),
                     annotation.value,
                     annotation.weight = 1)

  # return
  return(gxv)
}
