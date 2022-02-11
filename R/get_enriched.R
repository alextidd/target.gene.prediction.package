get_enriched <- function(variants,
                         H3K27ac,
                         specific_H3K27ac_closest_specific_genes,
                         contact,
                         expression,
                         TADs,
                         all_metadata,
                         out,
                         min_proportion_of_variants_in_top_H3K27ac,
                         tissue_of_interest,
                         do_all_cells,
                         do_all_cells_in_enriched_tissue,
                         estimate_cutoff = 2,
                         p.value_cutoff = 0.05){

  # for testing: # estimate_cutoff = 2 ; p.value_cutoff = 0.05

  enriched <- list()

  if(!is.null(tissue_of_interest)) {

    cat("Treating user-provided tissue of interest '", tissue_of_interest, "' as enriched tissue.\n")
    # User-provided tissue
    if(tissue_of_interest %ni% all_metadata$tissue){stop("Provided tissue of interest '", tissue_of_interest, "' not represented in the available data. Must be one of...\n",
                                                         paste(unique(all_metadata$tissue), collapse = ", "))}
    enriched[["celltypes"]] <- all_metadata %>%
      dplyr::filter(tissue == tissue_of_interest)

  } else if(do_all_cells) {

    cat("Applying annotations from all available cell types (`do_all_cells` = T).\n")
    # Include all cells
    enriched[["celltypes"]] <- all_metadata
    cat("All tissues: ") ; message(paste(unique(enriched$celltypes$tissue), collapse = ", "))
    cat("All celltypes: ") ; message(paste(unique(enriched$celltypes$name), collapse = ", "))

  } else {

    cat("Performing enrichment analysis to find enriched tissue(s).\n")
    # No user-provided tissue, determine via enrichment analysis
    # Fisher enrichment test will be of variants in top-decile cell-type-specificic H3K27ac marks in DHSs
    specific_H3K27ac <- H3K27ac$specificity %>%
      tidyr::pivot_longer(cols = -c(chrom:DHS),
                          names_to = "name",
                          values_to = "decile") %>%
      dplyr::filter(decile == 1)

    # threshold of % of CCVs in top H3K27ac
    counts <- bed_intersect_left(variants, specific_H3K27ac, keepBcoords = F) %>%
      dplyr::group_by(name) %>%
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
      bedA_groups = "name",
      bedB = variants,
      genome = ChrSizes) %>%
      dplyr::left_join(counts %>% dplyr::select(name, proportion_of_variants_in_top_H3K27ac), by = "name") %>%
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
      {dplyr::filter(all_metadata, name %in% .$name)} %>%
      # Get all samples in the enriched tissue, or only the enriched celltype
      {dplyr::filter(all_metadata,
                      (do_all_cells_in_enriched_tissue & tissue %in% .$tissue ) |
                      (!do_all_cells_in_enriched_tissue & name %in% .$name))}

    # Error message if no cell types were enriched
    if(nrow(enriched$celltypes)==0){stop("No enriched cell types found! Enrichment analysis saved to ", out$TissueEnrichment)}

    # Enriched celltype(s)/tissue(s)
    if(do_all_cells_in_enriched_tissue){
      cat("Enriched tissue(s): ") ; message(unique(enriched$celltypes$tissue))
      cat("Celltype(s) in enriched tissue(s): ") ; message(paste(unique(enriched$celltypes$name), collapse = ", "))
    } else {
      cat("Enriched celltype(s): ") ; message(paste(unique(enriched$celltypes$name), collapse = ", "))
    }

  }

  # Subset annotations
  ## H3K27ac
  enriched$H3K27ac <- H3K27ac %>%
    purrr::map( ~ dplyr::select(., chrom:DHS,
                                dplyr::any_of(enriched$celltypes$name[enriched$celltypes$object == "H3K27ac"])))
  ## specific_H3K27ac_closest_specific_genes
  enriched$specific_H3K27ac_closest_specific_genes <- specific_H3K27ac_closest_specific_genes %>%
    dplyr::select(chrom:end, dplyr::any_of(enriched$celltypes$name[enriched$celltypes$object == "H3K27ac"])) %>%
    dplyr::filter(dplyr::if_all(where(is.character), ~ !is.na(.)))
  ## contact
  enriched$contact <- names(contact) %>%
    sapply(function(x) {
      contact[[x]][names(contact[[x]]) %in% enriched$celltypes$name[enriched$celltypes$object == "contact"]]
    }, simplify = F, USE.NAMES = T)
  ## expression
  enriched$expression <- expression[, colnames(expression) %in% enriched$celltypes$name, drop = F]
  ## TADs
  enriched$TADs <- TADs[names(TADs) %in% enriched$celltypes$name]

  return(enriched)
}
