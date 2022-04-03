get_enriched <- function(variants,
                         DHSs,
                         H3K27ac_specificity_ranked,
                         H3K27ac,
                         expression,
                         expressed,
                         HiChIP,
                         TADs,
                         metadata,
                         out,
                         min_proportion_of_variants_in_top_H3K27ac,
                         celltype_of_interest,
                         tissue_of_interest,
                         do_all_celltypes,
                         do_all_celltypes_in_enriched_tissue,
                         ratio_cutoff = 1,
                         p_value_cutoff = 0.05){

  # for testing: # estimate_cutoff = 2 ; p_value_cutoff = 0.05

  enriched <- list()

  if (!is.null(tissue_of_interest)) {

    cat("Treating user-provided tissue of interest '", tissue_of_interest, "' as the enriched tissue.\n")
    # user-provided tissue
    enriched[["celltypes"]] <- metadata %>%
      dplyr::filter(tissue == tissue_of_interest)

  } else if (!is.null(celltype_of_interest)) {

    cat("Treating user-provided celltype of interest '", celltype_of_interest, "' as the enriched celltype.\n")
    # user-provided celltype
    enriched[["celltypes"]] <- metadata %>%
      dplyr::filter(celltype == celltype_of_interest)

  } else if (do_all_celltypes) {

    cat("Applying annotations from all available cell types (`do_all_celltypes` = T).\n")
    # Include all celltypes
    enriched[["celltypes"]] <- metadata
    cat("All tissues: ") ; enriched$celltypes$tissue %>% unique %>% paste(collapse = ", ") %>% cat('\n')
    cat("All celltypes: ") ; enriched$celltypes$celltype %>% unique %>% paste(collapse = ", ") %>% cat('\n')

  } else {

    cat("Performing enrichment analysis to find enriched tissue(s).\n")
    # No user-provided tissue, determine via enrichment analysis
    # Fisher enrichment test
    # Intersect SNPs with DHS sites and compute the mean H3K27ac specificity rank of the intersected DHS.
    # Assuming SNPs overlap sites randomly, compute significance of mean rank of sites where SNPs overlap based on deivation from uniform distribution.

    # intersect variants
    enrichment <- bed_intersect_left(variants, DHSs, keepBcoords = F) %>%
      {dplyr::filter(H3K27ac_specificity_ranked, DHS %in% .$DHS)} %>%
      # mean specificity rank per celltype
      dplyr::summarise(across(where(is.numeric), mean)) %>%
      tidyr::pivot_longer(everything(), names_to = "celltype", values_to = "mean_rank") %>%
      dplyr::mutate(
        # uniform distribution parameters
        N = nrow(DHSs),
        n = nrow(variants),
        mean = (N + 1)/2,
        variance = sqrt((N^2 - 1)/(12 * n)),
        # deviation from uniform distribution = enrichment
        ratio = mean_rank / mean,
        p_value = pnorm(mean_rank, mean, variance, lower.tail = F),
        p_value_adjust = p_value %>% p.adjust,
        pass = ((p_value_adjust < p_value_cutoff) & (ratio > ratio_cutoff))) %>%
      dplyr::arrange(p_value)
    write_tibble(enrichment, out$tissue_enrichment)

    # Enriched celltypes/tissues
    enriched[["celltypes"]] <- enrichment %>%
      # Filter to celltypes that pass filters
      dplyr::filter(pass)  %>%
      # Get enriched celltype metadata
      {dplyr::filter(metadata, celltype %in% .$celltype)} %>%
      # Get all samples in the enriched tissue, or only the enriched celltype
      {dplyr::filter(metadata,
                      (do_all_celltypes_in_enriched_tissue & tissue %in% .$tissue ) |
                      (!do_all_celltypes_in_enriched_tissue & celltype %in% .$celltype))}

    # Error message if no cell types were enriched
    if(nrow(enriched$celltypes)==0){
      stop("No enriched cell types found! Enrichment analysis saved to ", out$tissue_enrichment)}

    # Enriched celltype(s)/tissue(s)
    if(do_all_celltypes_in_enriched_tissue){
      cat("Enriched tissue(s): ") ; enriched$celltypes$tissue %>% unique %>% paste(collapse = ", ") %>% cat('\n')
      cat("Celltype(s) in enriched tissue(s): ") ; enriched$celltypes$celltype %>% unique %>% paste(collapse = ", ") %>% cat('\n')
    } else {
      cat("Enriched celltype(s): ") ; enriched$celltypes$celltype %>% unique %>% paste(collapse = ", ") %>% cat('\n')

      full_panel <- metadata$object
      enriched_panel <- dplyr::filter(metadata, celltype %in% enriched$celltypes$celltype)$object
      panel_gaps <- setdiff(full_panel, enriched_panel)
      if(length(panel_gaps) > 0){
      stop("\nOption `do_all_celltypes_in_enriched_tissue = F` was given, but the enriched celltype(s) do not constitute a full reference panel.",
           "\nEnriched celltype(s): ", enriched$celltypes$celltype  %>% unique %>% paste(collapse = ", "),
           "\nMissing annotations: ", paste(panel_gaps, collapse = ", "),
           "\nMake sure the enriched celltype has coverage across all annotations (TADs, HiChIP, expression, H3K27ac) in the metadata table.",
           "\nEither run again with `do_all_celltypes_in_enriched_tissue = T`` or lower the enrichment threshold `p_value_cutoff = ",p_value_cutoff,"`")
      }
    }

  }

  # Subset annotations
  ## H3K27ac
  enriched$H3K27ac <- H3K27ac %>%
    purrr::map( ~ dplyr::select(., DHS,
                                dplyr::any_of(enriched$celltypes$celltype[enriched$celltypes$object == "H3K27ac"])))
  ## HiChIP
  enriched$HiChIP <- HiChIP[names(HiChIP)[names(HiChIP) %in% enriched$celltypes$celltype[enriched$celltypes$object == "HiChIP"]]]
  ## expression
  enriched$expression <- expression %>%
    purrr::map( ~ dplyr::select(., ensg,
                                dplyr::any_of(enriched$celltypes$celltype[enriched$celltypes$object == "expression"])))
  ## expressed
  enriched$expressed <- expressed %>%
    dplyr::select(ensg, dplyr::any_of(enriched$celltypes$celltype[enriched$celltypes$object == "expression"]))
  ## TADs
  enriched$TADs <- TADs[names(TADs) %in% enriched$celltypes$celltype]

  return(enriched)
}


##### OLD ENRICHMENT TEST: #####
# # threshold of % of CCVs in top H3K27ac
# counts <- bed_intersect_left(variants, specific_H3K27ac, keepBcoords = F) %>%
#   dplyr::group_by(celltype) %>%
#   dplyr::count(name = "n_intersections") %>%
#   dplyr::mutate(n_variants = dplyr::n_distinct(variants$variant),
#                 proportion_of_variants_in_top_H3K27ac = n_intersections/n_variants)
# thresholded_counts <- counts %>%
#   dplyr::filter(proportion_of_variants_in_top_H3K27ac > min_proportion_of_variants_in_top_H3K27ac)
#
# # Fisher enrichment statistics for all celltypes/tissues
# enrichment <- bed_fisher_grouped(
#   bedA = specific_H3K27ac,
#   bedA_groups = "celltype",
#   bedB = variants,
#   genome = ChrSizes) %>%
#   dplyr::left_join(counts %>% dplyr::select(celltype, proportion_of_variants_in_top_H3K27ac), by = "celltype") %>%
#   dplyr::mutate(p_value.adjusted = p_value %>% p.adjust,
#                 pass =
#                   # filter for effect
#                   estimate > estimate_cutoff &
#                   # filter for significance
#                   p_value < p_value_cutoff &
#                   # threshold of % of CCVs in top H3K27ac
#                   proportion_of_variants_in_top_H3K27ac > min_proportion_of_variants_in_top_H3K27ac) %>%
#   dplyr::arrange(p_value)
# write_tibble(enrichment, out$tissue_enrichment)
#
# # Error message if no cell type exceeds min_proportion_of_variants_in_top_H3K27ac
# if(nrow(thresholded_counts) == 0){
#   max_prop <- counts %>%
#     dplyr::ungroup() %>%
#     dplyr::filter(proportion_of_variants_in_top_H3K27ac == max(proportion_of_variants_in_top_H3K27ac))
#   stop(message(
#     "No cell types had more than ",
#     min_proportion_of_variants_in_top_H3K27ac*100, "% of ", trait,
#     " variants in their most specific H3K27ac. Choose a lower `min_proportion_of_variants_in_top_H3K27ac` cut-off or specify a `tissue_of_interest` to skip this enrichment step.",
#     "\nCelltype(s) with the highest proportion = ", paste(max_prop$celltype, collapse = ", "),
#     "\nHighest proportion = ", unique(max_prop$proportion_of_variants_in_top_H3K27ac),
#     "\nEnrichment analysis saved to ", out$tissue_enrichment))
# }
