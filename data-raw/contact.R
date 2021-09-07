## code to prepare `contact` dataset goes here

contact <- list()

# Jung PC-HiC (ovary, lung)
# Must store repetitive metadata in list names, not in columns, to save space (cell type and assay name)
BEDPEFiles <- list.files("/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/PC-HiC/", full.names = T)
for(i in 1:length(BEDPEFiles)){
  ends <- target.gene.prediction.package::import_BEDPE_to_List(
    bedpefile = BEDPEFiles[i],
    metadata_cols = c("score", "CellType") # score = normalised_interaction_frequency
    )

  CellType <- unique(c(ends[[1]]$CellType, ends[[2]]$CellType)) ; cat("\n", CellType, "\n")
  annotation.name <-  paste0("PC-HiC_",CellType)

  ends <- purrr::map(ends, ~.x %>%
                       dplyr::select(-CellType))

  contact[[annotation.name]] <- ends
}

# Trench HiChIP (MCF7)
# Must store repetitive metadata in list names, not in columns, to save space (cell type and assay name)
contact[["HiChIP_breast"]] <- target.gene.prediction.package::import_BEDPE_to_List(
  bedpefile = "/working/lab_georgiat/alexandT/target_gene_prediction/data/HiChIP/MCF7/abc-ready/HiChIP.bedpe",
  metadata_cols = c("OldPairID", "score")
  ) %>%
  purrr::map(~.x %>%
               dplyr::select(-OldPairID) %>%
               # for infinite score values, set equal to the maximum non-infinite score
               dplyr::mutate(score = dplyr::case_when(is.infinite(score) ~ max(score[!is.infinite(score)]),
                                                      TRUE ~ score)) %>%
               # normalise to maximum
               dplyr::mutate(score = score / max(score)))

usethis::use_data(contact, overwrite = TRUE)

