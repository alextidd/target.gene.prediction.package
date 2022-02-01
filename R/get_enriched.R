get_enriched <- function(variants,
                         DHSs,
                         specific_DHSs_closest_specific_genes,
                         contact,
                         expression,
                         TADs,
                         all_metadata,
                         out,
                         min_proportion_of_variants_in_top_DHSs,
                         tissue_of_interest,
                         do_all_cells,
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
    specific_DHSs <- DHSs$specificity %>%
      tidyr::pivot_longer(cols = -c(chrom:DHS),
                          names_to = "name",
                          values_to = "decile") %>%
      dplyr::filter(decile == 1)

    # threshold of % of CCVs in top DHSs
    counts <- bed_intersect_left(variants, specific_DHSs, keepBcoords = F) %>%
      dplyr::group_by(name) %>%
      dplyr::count(name = "n_intersections") %>%
      dplyr::mutate(n_variants = dplyr::n_distinct(variants$variant),
                    proportion_of_variants_in_top_DHSs = n_intersections/n_variants)
    thresholded_counts <- counts %>%
      dplyr::filter(proportion_of_variants_in_top_DHSs > min_proportion_of_variants_in_top_DHSs)
    if(nrow(thresholded_counts) == 0){
      stop(message("No cell types had more than ",
                   min_proportion_of_variants_in_top_DHSs*100, "% of ", trait,
                   " variants in their most specific DHSs. Choose a lower `min_proportion_of_variants_in_top_DHSs` cut-off or specify a `tissue_of_interest` to skip this enrichment step."))
    }

    # Fisher enrichment statistics for all celltypes/tissues
    enrichment <- bed_fisher_grouped(
      bedA = specific_DHSs,
      bedA_groups = "name",
      bedB = variants,
      genome = ChrSizes) %>%
      dplyr::left_join(counts %>% dplyr::select(name, proportion_of_variants_in_top_DHSs)) %>%
      dplyr::mutate(enriched =
                      # filter for effect
                      estimate > estimate_cutoff &
                      # filter for significance
                      p.value < p.value_cutoff &
                      # threshold of % of CCVs in top DHSs
                      proportion_of_variants_in_top_DHSs > min_proportion_of_variants_in_top_DHSs)
    write_tibble(enrichment, out$TissueEnrichment)

    # Enriched celltypes/tissues
    enriched[["celltypes"]] <- enrichment %>%
      # Filter to enriched celltypes
      dplyr::filter(enriched)  %>%
      # Get enriched celltype metadata
      {dplyr::filter(all_metadata, name %in% .)} %>%
      # Get all samples in the enriched tissue
      {dplyr::filter(all_metadata, tissue %in% .$tissue)}

    # Error message if no cell types were enriched
    if(nrow(enriched$celltypes)==0){stop("No enriched cell types found!")}

    # Enriched celltype(s)/tissue(s)
    cat("Enriched tissue(s): ") ; message(unique(enriched$celltypes$tissue))
    cat("Celltype(s) in enriched tissue(s): ") ; message(paste(unique(enriched$celltypes$name), collapse = ", "))

  }

  # Subset annotations
  ## DHSs
  enriched$DHSs <- DHSs %>%
    purrr::map( ~ dplyr::select(., chrom:DHS,
                                dplyr::any_of(enriched$celltypes$name[enriched$celltypes$object == "DHSs"])))
  ## specific_DHSs_closest_specific_genes
  enriched$specific_DHSs_closest_specific_genes <-
    specific_DHSs_closest_specific_genes %>%
    dplyr::select(chrom:end, dplyr::all_of(enriched$celltypes$name[enriched$celltypes$object == "DHSs"])) %>%
    dplyr::filter(dplyr::if_all(where(is.character), ~ !is.na(.)))
  ## contact
  enriched$contact <-
    names(contact) %>% sapply(function(x) {
      contact[[x]][names(contact[[x]]) %in% enriched$celltypes$name[enriched$celltypes$object == "contact"]]
    },
    simplify = F, USE.NAMES = T)
  # enriched$contact <-
  #   enriched$contact[lapply(enriched$contact, length) > 0] # remove empty elements
  ## expression
  enriched$expression <-
    expression[, colnames(expression) %in% enriched$celltypes$name, drop = F]
  ## TADs
  enriched$TADs <- TADs[names(TADs) %in% enriched$celltypes$name]

  return(enriched)
}
