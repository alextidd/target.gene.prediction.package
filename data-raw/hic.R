## code to prepare `hichip` dataset goes here

BEDPEFiles <- list.files("/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/PC-HiC/", full.names = T)

hic <- list()
for(i in 1:length(BEDPEFiles)){
  ends <- target.gene.prediction.package::import_BEDPE_to_List(
    bedpefile = BEDPEFiles[i],
    metadata_cols = c("CellType", "normalised_interaction_frequency"))

  ends <- purrr::map(ends, ~dplyr::filter(.x, normalised_interaction_frequency > 0))

  CellType <- unique(c(ends[[1]]$CellType, ends[[2]]$CellType)) ; cat(CellType, "\n")
  hic[[CellType]] <- ends
}

usethis::use_data(hic, overwrite = TRUE)


# # HiChIP example dataset
# BEDPEFile <- "/working/lab_georgiat/alexandT/target_gene_prediction/data/HiChIP/MCF7/abc-ready/HiChIP.bedpe"
# hichip <- target.gene.prediction.package::import_BEDPE_to_List(bedpefile = BEDPEFile,
#                                                                metadata_cols = c("OldPairID", "value"))
