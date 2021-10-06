get_txv_level_annotations <- function(open_variants,
                                      variant_to_gene_max_distance,
                                      enriched,
                                      contact,
                                      TADs) {
  cat("Annotating transcript x variant pairs...\n")

  txv <- list()

  # TSS distance scoring method (these pairings are the only ones to consider)
  distance <- open_variants %>%
    # Get genes within variant_to_gene_max_distance of each variant
    valr::bed_slop(both = variant_to_gene_max_distance,
                   genome = target.gene.prediction.package::ChrSizes,
                   trim = T) %>%
    valr::bed_intersect(., target.gene.prediction.package::TSSs,
                        suffix = c(".variant", ".TSS")) %>%
    dplyr::group_by(variant.variant) %>%
    # Calculate inverse of the absolute bp distance for each variant-transcript pair
    dplyr::mutate(pair.distance = abs((start.variant + variant_to_gene_max_distance) - start.TSS),
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
  txv$inv_distance <- distance %>%
    dplyr::transmute(variant, enst,
                     value = pair.inverse_distance)

  # variant-TSS distance rank method
  txv$inv_distance_rank <- distance %>%
    dplyr::transmute(variant, enst,
                     value = pair.inverse_distance_rank)

  # closest variant-TSS method
  txv$closest <- distance %>%
    dplyr::filter(pair.inverse_distance_rank == 1) %>%
    dplyr::transmute(variant, enst,
                     value = 1)

  # intersect loop ends, by cell type, with enhancer variants and gene TSSs
  # (finds interaction loops with a variant at one end and a TSS at the other)
  txv <- contact %>%
    # Intersect with the contact data
    purrr::map(~ intersect_BEDPE(
      # ! For mutually exclusive intersection with HiC ranges, make variant intervals 1bp long, equal to the end position
      SNPend = open_variants %>% dplyr::mutate(start = end),
      TSSend = target.gene.prediction.package::TSSs,
      bedpe = .) %>%
        tidyr::separate(InteractionID, into = c("celltype", "assay"), sep = "\\_", remove = F, extra = "drop") %>%
        dplyr::transmute(variant, enst,
                         InteractionID,
                         value = score,
                         celltype,
                         assay = paste0("contact_", assay))
      ) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    # Make sure all interactions are within 2Mb - hard filter, ignore everything further
    dplyr::inner_join(distance %>% dplyr::select(variant, enst)) %>%
    # Split into different contact assays
    split(f = .$assay) %>%
    # Widen
    purrr::map(~ tidyr::pivot_wider(.,
      id_cols = c(variant, enst),
      names_from = celltype,
      values_from = value)) %>%
    # Add to txv list
    c(., txv)

  ## ~10% of variant-TSS interactions indicated by the contact data
  ## are further than 2Mb apart and are thus eliminated

  # get variants at promoters - score by sum of signal and specificity per enriched celltype at the DHS (~promoter activity)
  txv$promoter <- open_variants %>%
    dplyr::select(chrom:variant) %>%
    # Get variants within promoter regions
    target.gene.prediction.package::bed_intersect_left(
      target.gene.prediction.package::promoters,
      keepBcoords = F, keepBmetadata = T) %>%
    # Get DHS bins at those promoter variants
    target.gene.prediction.package::intersect_DHSs(list(), ., enriched$DHSs) %>%
    dplyr::bind_rows() %>%
    # Sum specificity + signal bin per promoter variant per enriched cell type
    dplyr::group_by(dplyr::across(setdiff(names(.), enriched$tissues$code))) %>%
    dplyr::summarise(dplyr::across(everything(), sum)) %>%
    dplyr::ungroup()

  # intronic variants
  txv$intron <- open_variants %>%
    # Get variants within introns
    target.gene.prediction.package::bed_intersect_left(
      target.gene.prediction.package::introns, .,
      keepBcoords = F, keepBmetadata = T) %>%
    dplyr::transmute(variant,
                     enst,
                     value = 1)

  # exonic (coding) variants
  txv$exon <- open_variants %>%
    # Get variants within introns
    target.gene.prediction.package::bed_intersect_left(
      target.gene.prediction.package::exons, .,
      keepBcoords = F, keepBmetadata = T) %>%
    dplyr::transmute(variant,
                     enst,
                     value = 1)

  # TADs
  txv$TADs <- dplyr::full_join(
    target.gene.prediction.package::bed_intersect_left(
      open_variants, TADs, keepBcoords = F) %>%
      dplyr::select(variant, TAD),
    target.gene.prediction.package::bed_intersect_left(
      target.gene.prediction.package::TSSs, TADs, keepBcoords = F) %>%
      dplyr::select(enst, TAD)
  ) %>%
    dplyr::transmute(variant, enst, value = 1)

  # return
  names(txv) <- paste0("txv_", names(txv))
  return(txv)
}
