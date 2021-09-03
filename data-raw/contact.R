## code to prepare `hichip` dataset goes here

contact <- list()

BEDPEFiles <- list.files("/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/PC-HiC/", full.names = T)
for(i in 1:length(BEDPEFiles)){
  ends <- target.gene.prediction.package::import_BEDPE_to_List(
    bedpefile = BEDPEFiles[i],
    metadata_cols = c("CellType", "score") # score = normalised_interaction_frequency
    )

  CellType <- unique(c(ends[[1]]$CellType, ends[[2]]$CellType)) ; cat(CellType, "\n")
  annotation.name <-  paste0("PC-HiC_",CellType)

  ends <- purrr::map(ends, ~.x %>%
                       dplyr::filter(score > 0) %>%
                       dplyr::mutate(annotation.name = annotation.name))

  contact[[annotation.name]] <- ends
}

# MCF7
contact[["HiChIP_breast"]] <- target.gene.prediction.package::import_BEDPE_to_List(
  bedpefile = "/working/lab_georgiat/alexandT/target_gene_prediction/data/HiChIP/MCF7/abc-ready/HiChIP.bedpe",
  metadata_cols = c("OldPairID", "score")
  ) %>%
  purrr::map(~.x %>%
               dplyr::select(-OldPairID) %>%
               dplyr::mutate(CellType = "breast",
                             annotation.name = "HiChIP_breast"))

usethis::use_data(contact, overwrite = TRUE)




# # HiChIP example dataset
# BEDPEFile <- "/working/lab_georgiat/alexandT/target_gene_prediction/data/HiChIP/MCF7/abc-ready/HiChIP.bedpe"
# hichip <- target.gene.prediction.package::import_BEDPE_to_List(bedpefile = BEDPEFile,
#                                                                metadata_cols = c("OldPairID", "value"))
