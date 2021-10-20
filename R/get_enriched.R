get_enriched <- function(DHSs,
                         DHSs_metadata,
                         contact_metadata,
                         variants,
                         min_proportion_of_variants_in_top_DHSs,
                         tissue_of_interest,
                         estimate_cutoff = 2,
                         p.value_cutoff = 0.05){

  # for testing:
  # estimate_cutoff = 2 ; p.value_cutoff = 0.05

  enriched <- list()

  # metadata for all annotations
  all_metadata <- DHSs_metadata %>% dplyr::mutate(object = "DHSs") %>%
    dplyr::bind_rows(contact_metadata %>% dplyr::mutate(object = "contact"))

  if(!is.null(tissue)){
    # User-provided tissue
    if(tissue %ni% all_metadata$tissue){stop("Provided tissue '", tissue, "' not represented in the available data. Must be one of...\n",
                                            paste(unique(all_metadata$tissue), collapse = ", "))}

    cat("Treating user-provided tissue '", tissue, "' as enriched tissue.\n")
    enriched[["tissues"]] <- all_metadata %>%
      dplyr::filter(tissue == tissue_of_interest)

  } else {
    # No user-provided tissue, find via enrichment analysis

    # Fisher enrichment test will be of variants in upper-quartile cell-type-specificic H3K27ac marks in DHSs
    specific_DHSs <- DHSs$specificity %>%
      tidyr::gather(key = "annotation",
                    value = "decile",
                    -c(chrom:DHS)) %>%
      dplyr::filter(decile == 1)

    # threshold of % of CCVs in top DHSs
    thresholded_counts <- target.gene.prediction.package::bed_intersect_left(variants, specific_DHSs, keepBcoords = F) %>%
      dplyr::group_by(annotation) %>%
      dplyr::count(name = "n_intersections") %>%
      dplyr::mutate(n_variants = dplyr::n_distinct(variants$variant)) %>%
      dplyr::filter(n_intersections/n_variants > min_proportion_of_variants_in_top_DHSs)

    enriched[["celltypes"]] <- target.gene.prediction.package::bed_fisher_grouped(
      bedA = specific_DHSs,
      bedA_groups = "annotation",
      bedB = variants,
      genome = target.gene.prediction.package::ChrSizes,
      # filter for effect and significance
      estimate > estimate_cutoff,
      p.value < p.value_cutoff
    ) %>%
      # threshold of % of CCVs in top DHSs
      dplyr::filter(annotation %in% thresholded_counts$annotation) %>%
      # Extract enriched cell types
      dplyr::pull(annotation) %>%
      {dplyr::filter(DHSs_metadata, mnemonic %in% .)}
    enriched[["tissues"]] <- all_metadata %>%
      dplyr::filter(tissue %in% enriched$celltypes$tissue)

    cat("Enriched cell type(s): ", enriched$celltypes$mnemonic, "\n")
  }
  cat("Enriched tissue(s):", unique(enriched$tissues$tissue), "\n")

  return(enriched)
}
