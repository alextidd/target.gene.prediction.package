get_enriched <- function(variants,
                         H3K27ac,
                         specific_H3K27ac_closest_specific_genes,
                         contact,
                         expression,
                         expressed,
                         TADs,
                         all_metadata,
                         out,
                         min_proportion_of_variants_in_top_H3K27ac,
                         celltype_of_interest,
                         tissue_of_interest,
                         do_all_celltypes,
                         do_all_celltypes_in_enriched_tissue,
                         estimate_cutoff = 2,
                         p.value_cutoff = 0.05){

  # for testing: # estimate_cutoff = 2 ; p.value_cutoff = 0.05

  enriched <- list()

  if (!is.null(tissue_of_interest)) {

    cat("Treating user-provided tissue of interest '", tissue_of_interest, "' as the enriched tissue.\n")
    # user-provided tissue
    enriched[["celltypes"]] <- all_metadata %>%
      dplyr::filter(tissue == tissue_of_interest)

  } else if (!is.null(celltype_of_interest)) {

    cat("Treating user-provided celltype of interest '", celltype_of_interest, "' as the enriched celltype.\n")
    # user-provided celltype
    enriched[["celltypes"]] <- all_metadata %>%
      dplyr::filter(celltype == celltype_of_interest)

  } else if (do_all_celltypes) {

    cat("Applying annotations from all available cell types (`do_all_celltypes` = T).\n")
    # Include all celltypes
    enriched[["celltypes"]] <- all_metadata
    cat("All tissues: ") ; message(paste(unique(enriched$celltypes$tissue), collapse = ", "))
    cat("All celltypes: ") ; message(paste(unique(enriched$celltypes$celltype), collapse = ", "))

  } else {

    cat("Performing enrichment analysis to find enriched tissue(s).\n")
    # No user-provided tissue, determine via enrichment analysis
    # Fisher enrichment test will be of variants in top-decile cell-type-specificic H3K27ac marks in DHSs
    specific_H3K27ac <- H3K27ac$specificity %>%
      tidyr::pivot_longer(cols = -c(chrom:DHS),
                          names_to = "celltype",
                          values_to = "decile") %>%
      dplyr::filter(decile == 1)

    # threshold of % of CCVs in top H3K27ac
    counts <- bed_intersect_left(variants, specific_H3K27ac, keepBcoords = F) %>%
      dplyr::group_by(celltype) %>%
      dplyr::count(name = "n_intersections") %>%
      dplyr::mutate(n_variants = dplyr::n_distinct(variants$variant),
                    proportion_of_variants_in_top_H3K27ac = n_intersections/n_variants)
    thresholded_counts <- counts %>%
      dplyr::filter(proportion_of_variants_in_top_H3K27ac > min_proportion_of_variants_in_top_H3K27ac)
    if(nrow(thresholded_counts) == 0){
      stop(message("No cell types had more than ",
                   min_proportion_of_variants_in_top_H3K27ac*100, "% of ", trait,
                   " variants in their most specific H3K27ac. Choose a lower `min_proportion_of_variants_in_top_H3K27ac` cut-off or specify a `tissue_of_interest` to skip this enrichment step."))
    }

    # Fisher enrichment statistics for all celltypes/tissues
    enrichment <- bed_fisher_grouped(
      bedA = specific_H3K27ac,
      bedA_groups = "celltype",
      bedB = variants,
      genome = ChrSizes) %>%
      dplyr::left_join(counts %>% dplyr::select(celltype, proportion_of_variants_in_top_H3K27ac), by = "celltype") %>%
      dplyr::mutate(significant =
                      # filter for effect
                      estimate > estimate_cutoff &
                      # filter for significance
                      p.value < p.value_cutoff &
                      # threshold of % of CCVs in top H3K27ac
                      proportion_of_variants_in_top_H3K27ac > min_proportion_of_variants_in_top_H3K27ac) %>%
      dplyr::arrange(p.value)
    write_tibble(enrichment, out$TissueEnrichment)

    # Enriched celltypes/tissues
    enriched[["celltypes"]] <- enrichment %>%
      # Filter to enriched celltypes
      dplyr::filter(significant)  %>%
      # Get enriched celltype metadata
      {dplyr::filter(all_metadata, celltype %in% .$celltype)} %>%
      # Get all samples in the enriched tissue, or only the enriched celltype
      {dplyr::filter(all_metadata,
                      (do_all_celltypes_in_enriched_tissue & tissue %in% .$tissue ) |
                      (!do_all_celltypes_in_enriched_tissue & celltype %in% .$celltype))}

    # Error message if no cell types were enriched
    if(nrow(enriched$celltypes)==0){stop("No enriched cell types found! Enrichment analysis saved to ", out$TissueEnrichment)}

    # Enriched celltype(s)/tissue(s)
    if(do_all_celltypes_in_enriched_tissue){
      cat("Enriched tissue(s): ") ; enriched$celltypes$tissue %>% unique %>% paste(collapse = ", ") %>% cat('\n')
      cat("Celltype(s) in enriched tissue(s): ") ; enriched$celltypes$celltype %>% unique %>% paste(collapse = ", ") %>% cat('\n')
    } else {
      cat("Enriched celltype(s): ") ; enriched$celltypes$celltype %>% unique %>% paste(collapse = ", ") %>% cat('\n')
    }

  }

  # Subset annotations
  ## H3K27ac
  enriched$H3K27ac <- H3K27ac %>%
    purrr::map( ~ dplyr::select(., chrom:DHS,
                                dplyr::any_of(enriched$celltypes$celltype[enriched$celltypes$object == "H3K27ac"])))
  ## specific_H3K27ac_closest_specific_genes
  enriched$specific_H3K27ac_closest_specific_genes <- specific_H3K27ac_closest_specific_genes %>%
    dplyr::select(chrom:end, dplyr::any_of(enriched$celltypes$celltype[enriched$celltypes$object == "H3K27ac"])) %>%
    dplyr::filter(dplyr::if_all(where(is.character), ~ !is.na(.)))
  ## contact
  enriched$contact <- names(contact) %>%
    sapply(function(x) {
      contact[[x]][names(contact[[x]]) %in% enriched$celltypes$celltype[enriched$celltypes$object == "contact"]]
    }, simplify = F, USE.NAMES = T)
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
