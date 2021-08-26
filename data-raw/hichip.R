## code to prepare `hichip` dataset goes here

# HiChIP example dataset
BEDPEFile <- "/working/lab_georgiat/alexandT/target_gene_prediction/data/HiChIP/MCF7/abc-ready/HiChIP.bedpe"
hichip <- target.gene.prediction.package::import_BEDPE_to_List(bedpefile = BEDPEFile,
                                                               metadata_cols = c("OldPairID", "value"))

usethis::use_data(hichip, overwrite = TRUE)
