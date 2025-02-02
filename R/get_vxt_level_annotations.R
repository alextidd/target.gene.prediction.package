get_vxt_level_annotations <- function(variants,
                                      DHSs,
                                      vxt_master,
                                      variant_to_gene_max_distance,
                                      enriched) {

  vxt <- list()

  # TSS distance scoring method (these pairings are the only ones to consider)
  distance <- vxt_master %>%
    dplyr::group_by(variant) %>%
    # Calculate inverse of the absolute bp distance for each variant-transcript pair (avoid infinite values)
    dplyr::mutate(inv_distance = dplyr::case_when(distance == 0 ~ 1,
                                                  TRUE ~ 1 / distance),
                  # ranking transcript TSSs (if two transcript TSSs are equidistant to the variant, they will receive the same, lower rank)
                  inv_distance_rank = 1 / rank(distance, ties.method = "min"))

  # variant-TSS inverse distance score method
  vxt$inv_distance <- distance %>%
    dplyr::transmute(cs, variant, enst,
                     value = inv_distance)

  # variant-TSS distance rank method
  vxt$inv_distance_rank <- distance %>%
    dplyr::transmute(cs, variant, enst,
                     value = inv_distance_rank)

  # closest variant-TSS method
  vxt$closest <- distance %>%
    dplyr::filter(inv_distance_rank == 1) %>%
    dplyr::transmute(cs, variant, enst)

  # intersect loop ends, by cell type, with enhancer variants and gene TSSs
  # (finds interaction loops with a variant at one end and a TSS at the other)
  vxt$HiChIP_scores <- enriched$HiChIP %>% names %>%
    sapply(function(celltype){
        # Intersect with the HiChIP data
        enriched$HiChIP[[celltype]] %>%
          intersect_BEDPE(
            # ! For mutually exclusive intersection with ranges, make variant intervals 1bp long, equal to the end position
            SNPend = variants %>% dplyr::mutate(start = end),
            TSSend = TSSs,
            bedpe = .) %>%
          dplyr::transmute(cs, variant, enst,
                           InteractionID,
                           value = score) %>%
          dplyr::inner_join(distance %>% dplyr::select(cs, variant, enst), by = c("cs", "variant", "enst"))
      }, simplify = F, USE.NAMES = T) %>%
        dplyr::bind_rows(.id = "celltype") %>%
        #dplyr::group_by(celltype, variant, enst) %>%
        #dplyr::count()
        #dplyr::summarise(value = sum(value)) %>%
        tidyr::pivot_wider(id_cols = c(cs, variant, enst),
                           names_from = celltype,
                           values_from = value)

  # HiChIP binary
  vxt$HiChIP_binary <- vxt$HiChIP_scores %>%
    dplyr::mutate(., dplyr::across(where(is.numeric), ~ dplyr::case_when(is.na(.) ~ 0, TRUE ~ 1)))

  # get variants at promoters
  vxt$promoter <- variants %>%
    dplyr::select(chrom:end, variant, cs) %>%
    # Get variants within promoter regions
    bed_intersect_left(
      promoters,
      keepBcoords = F, keepBmetadata = T) %>%
    dplyr::transmute(cs, variant, enst)

  # get variants at promoters - score by sum of signal and specificity per celltype at the DHS (~promoter activity)
  vxt$promoter_H3K27ac_bins_sum <- variants %>%
    dplyr::select(chrom:end, variant, cs) %>%
    # Get variants within promoter regions
    bed_intersect_left(
      promoters,
      keepBcoords = F, keepBmetadata = T) %>%
    # Get DHS bins at those promoter variants
    intersect_H3K27ac(list(),
                   query = .,
                   DHSs,
                   H3K27ac = enriched$H3K27ac,
                   cs, variant, enst) %>%
    purrr::reduce(dplyr::bind_rows) %>%
    # Sum specificity + signal bin per promoter variant per cell type
    dplyr::group_by(cs, variant, enst) %>%
    dplyr::summarise(dplyr::across(where(is.numeric), sum)) %>%
    dplyr::ungroup()

  # intronic variants
  vxt$intron <- variants %>%
    # Get variants within introns
    bed_intersect_left(
      introns, .,
      keepBcoords = F, keepBmetadata = T) %>%
    dplyr::transmute(cs, variant, enst)

  # exonic (coding) variants
  vxt$exon <- variants %>%
    # Get variants within exons
    bed_intersect_left(
      exons, .,
      keepBcoords = F, keepBmetadata = T) %>%
    dplyr::transmute(cs, variant, enst)

  # TADs
  TADs_w_ID <- enriched$TADs %>%
    dplyr::bind_rows(.id = "celltype") %>%
    dplyr::group_by(celltype) %>%
    dplyr::mutate(TAD = paste0(celltype, "_", dplyr::row_number())) %>%
    dplyr::ungroup()
  vxt$TADs <- dplyr::inner_join(
    bed_intersect_left(
      variants, TADs_w_ID, keepBcoords = F) %>%
      dplyr::select(cs, variant, TAD, celltype),
    bed_intersect_left(
      TSSs, TADs_w_ID, keepBcoords = F) %>%
      dplyr::select(enst, TAD, celltype),
    by = c("TAD", "celltype")
  ) %>%
    dplyr::transmute(cs, variant, enst, celltype, value = 1) %>%
    tidyr::pivot_wider(id_cols = c("cs", "variant", "enst"),
                       names_from = celltype,
                       values_from = value)

  # return
  names(vxt) <- paste0("vxt_", names(vxt))
  return(vxt)
}
