## code to prepare `TSSs` dataset goes here

# Import a BED file of transcription start sites

TSSFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/reference_data/GENCODE/proteincoding.TSSs.GENCODE.bed.gz"
TSSs <- target.gene.prediction.package::import_BED(gzfile(TSSFile),
                                                   metadata_cols = c("ensg", "symbol", "enst"))

usethis::use_data(TSSs, overwrite = TRUE)
