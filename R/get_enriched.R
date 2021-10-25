get_enriched <- function(variants,
                         DHSs,
                         specific_DHSs_closest_specific_genes,
                         contact,
                         expression,
                         all_metadata,
                         min_proportion_of_variants_in_top_DHSs,
                         tissue_of_interest,
                         do_all_cells,
                         estimate_cutoff = 2,
                         p.value_cutoff = 0.05){

  # for testing:
  # estimate_cutoff = 2 ; p.value_cutoff = 0.05

  enriched <- list()

  if(!is.null(tissue_of_interest)){
    # User-provided tissue
    if(tissue_of_interest %ni% all_metadata$tissue){stop("Provided tissue of interest '", tissue_of_interest, "' not represented in the available data. Must be one of...\n",
                                            paste(unique(all_metadata$tissue), collapse = ", "))}

    cat("Treating user-provided tissue of interest '", tissue_of_interest, "' as enriched tissue.\n")
    enriched[["tissues"]] <- all_metadata %>%
      dplyr::filter(tissue == tissue_of_interest)

  } else {
    # No user-provided tissue, determine via enrichment analysis

    # Fisher enrichment test will be of variants in upper-quartile cell-type-specificic H3K27ac marks in DHSs
    specific_DHSs <- DHSs$specificity %>%
      tidyr::gather(key = "name",
                    value = "decile",
                    -c(chrom:DHS)) %>%
      dplyr::filter(decile == 1)

    # threshold of % of CCVs in top DHSs
    thresholded_counts <- target.gene.prediction.package::bed_intersect_left(variants, specific_DHSs, keepBcoords = F) %>%
      dplyr::group_by(name) %>%
      dplyr::count(name = "n_intersections") %>%
      dplyr::mutate(n_variants = dplyr::n_distinct(variants$variant)) %>%
      dplyr::filter(n_intersections/n_variants > min_proportion_of_variants_in_top_DHSs)

    enriched[["celltypes"]] <- target.gene.prediction.package::bed_fisher_grouped(
      bedA = specific_DHSs,
      bedA_groups = "name",
      bedB = variants,
      genome = target.gene.prediction.package::ChrSizes,
      # filter for effect and significance
      estimate > estimate_cutoff,
      p.value < p.value_cutoff
    ) %>%
      # threshold of % of CCVs in top DHSs
      dplyr::filter(name %in% thresholded_counts$name) %>%
      # Extract enriched cell types
      dplyr::pull(name) %>%
      {dplyr::filter(all_metadata, name == ., object == "DHSs")}
    enriched[["tissues"]] <- all_metadata %>%
      dplyr::filter(tissue %in% enriched$celltypes$tissue)
    # Error message if no cell types were enriched
    if(nrow(enriched$celltypes)==0){stop("No enriched cell types found!")}

    cat("Enriched cell type(s): ", enriched$celltypes$name, "\n")
  }
  cat("Enriched tissue(s):", unique(enriched$tissues$tissue), "\n")


  # Subset annotations to those in enriched tissues, or consider all cell types (do_all_cells argument, default = F)
  if(do_all_cells){
    cat("Applying annotations in all cell types...\n")
    enriched$DHSs <- DHSs
    enriched$specific_DHSs_closest_specific_genes <- specific_DHSs_closest_specific_genes
    enriched$contact <- contact
    enriched$expression <- expression
  } else {
    cat("Applying annotations in enriched cell type(s) only...\n")
    enriched$DHSs <- DHSs %>%
      purrr::map(~ dplyr::select(., chrom:DHS, dplyr::any_of(enriched$tissues$name[enriched$tissues$object == "DHSs"])))
    enriched$specific_DHSs_closest_specific_genes <- specific_DHSs_closest_specific_genes %>%
      dplyr::select(DHS, dplyr::any_of(enriched$tissues$name[enriched$tissues$object == "DHSs"]))
    enriched$contact <- names(contact) %>% sapply(function(x) {contact[[x]][names(contact[[x]]) %in% enriched$tissues$name[enriched$tissues$object == "contact"]]},
                                                  simplify = F, USE.NAMES = T)
    enriched$contact <- enriched$contact[lapply(enriched$contact, length) > 0] # remove empty elements
    enriched$expression <- expression[, colnames(expression) %in% enriched$tissues$name, drop = F]
  }

  return(enriched)
}
