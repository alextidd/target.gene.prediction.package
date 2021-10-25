get_txv_level_annotations <- function(variants,
                                      txv_master,
                                      variant_to_gene_max_distance,
                                      enriched,
                                      TADs) {
  cat("Annotating transcript x variant pairs...\n")

  txv <- list()

  # TSS distance scoring method (these pairings are the only ones to consider)
  distance <- txv_master %>%
    dplyr::group_by(variant) %>%
    # Calculate inverse of the absolute bp distance for each variant-transcript pair
    dplyr::mutate(pair.distance = abs((start.variant + variant_to_gene_max_distance) - start.TSS),
                  pair.inverse_distance = 1/pair.distance,
                  # ranking transcript TSSs (if two transcript TSSs are equidistant to the variant, they will receive the same, lower rank)
                  pair.inverse_distance_rank = 1/rank(pair.distance, ties.method = "min"))

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
  txv_contact_scores <- enriched$contact %>%
    purrr::map(~
    # Intersect with the contact data
    purrr::map(., ~ intersect_BEDPE(
      # ! For mutually exclusive intersection with ranges, make variant intervals 1bp long, equal to the end position
      SNPend = variants %>% dplyr::mutate(start = end),
      TSSend = target.gene.prediction.package::TSSs,
      bedpe = .) %>%
        dplyr::transmute(variant, enst,
                         InteractionID,
                         value = score) %>%
        dplyr::inner_join(distance %>% dplyr::select(variant, enst))
      ) %>%
      dplyr::bind_rows(.id = "celltype") %>%
      tidyr::pivot_wider(id_cols = c(variant, enst),
                         names_from = celltype,
                         values_from = value))
  names(txv_contact_scores) <- paste0("contact_", names(txv_contact_scores))
  txv <- c(txv, txv_contact_scores)

  # contact binary
  txv_contact_binary <- txv_contact_scores %>%
    purrr::map(~ dplyr::mutate(., dplyr::across(where(is.numeric), ~ dplyr::case_when(is.na(.) ~ 0,
                                                                                      TRUE ~ 1))))
  names(txv_contact_binary) <- paste0(names(txv_contact_binary), "_", "binary")
  txv <- c(txv, txv_contact_binary)

  # Problem to FIX!!! 4 datasets in which the bulk (>70%) of the distribution of scores is equal to the minimum score
  # - colorectal_HiChIP
  # - Hct116_ChIAPET
  # - Helas3_ChIAPET
  # - Nb4_ChIAPET
  # contact %>% names %>%
  #   purrr::map(~ data.frame(assay = .,
  #                           n.scores = nrow(contact[[.]]$first),
  #                           n.above.or.equal.median.score = dplyr::filter(contact[[.]]$first, score >= median(score)) %>% nrow,
  #                           n.above.median.score = dplyr::filter(contact[[.]]$first, score > median(score)) %>% nrow,
  #                           n.equal.median.score = dplyr::filter(contact[[.]]$first, score == median(score)) %>% nrow,
  #                           n.equal.min.score = dplyr::filter(contact[[.]]$first, score == min(score)) %>% nrow,
  #                           n.distinct.scores = contact[[.]]$first$score %>% dplyr::n_distinct()))  %>%
  #   purrr::reduce(dplyr::bind_rows) %>% tibble::as_tibble() %>%
  #   dplyr::mutate(percent.of.scores.equal.to.median = (n.equal.median.score/n.scores)*100,
  #                 percent.of.scores.equal.to.min = (n.equal.min.score/n.scores)*100) %>%
  #   dplyr::filter(n.above.or.equal.median.score > ((n.scores/2) + (0.01*n.scores)))

  ## ~10% of variant-TSS interactions indicated by the contact data
  ## are further than 2Mb apart and are thus eliminated

  # get variants at promoters
  txv$promoter <- variants %>%
    dplyr::select(chrom:variant) %>%
    # Get variants within promoter regions
    target.gene.prediction.package::bed_intersect_left(
      target.gene.prediction.package::promoters,
      keepBcoords = F, keepBmetadata = T) %>%
    dplyr::transmute(variant,
                     enst,
                     value = 1)

  # get variants at promoters - score by sum of signal and specificity percelltype at the DHS (~promoter activity)
  txv$promoter_DHS_bins_sum <- variants %>%
    dplyr::select(chrom:variant) %>%
    # Get variants within promoter regions
    target.gene.prediction.package::bed_intersect_left(
      target.gene.prediction.package::promoters,
      keepBcoords = F, keepBmetadata = T) %>%
    # Get DHS bins at those promoter variants
    intersect_DHSs(list(),
                   query = .,
                   DHSs = enriched$DHSs,
                   variant, enst) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    # Sum specificity + signal bin per promoter variant per cell type
    dplyr::group_by(variant, enst) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), sum)) %>%
    dplyr::ungroup()

  # intronic variants
  txv$intron <- variants %>%
    # Get variants within introns
    target.gene.prediction.package::bed_intersect_left(
      target.gene.prediction.package::introns, .,
      keepBcoords = F, keepBmetadata = T) %>%
    dplyr::transmute(variant,
                     enst,
                     value = 1)

  # exonic (coding) variants
  txv$exon <- variants %>%
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
      variants, TADs, keepBcoords = F) %>%
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
