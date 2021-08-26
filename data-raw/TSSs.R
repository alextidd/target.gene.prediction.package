## code to prepare `TSSs` dataset goes here

# Import a BED file of transcription start sites

TSSFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/reference_data/GENCODE/gencode.v34lift37.basic.tss.bed.gz"
TSSs <- target.gene.prediction.package::import_BED(gzfile(TSSFile),
                                                   metadata_cols = c("ensg", "symbol", "enst"))

usethis::use_data(TSSs, overwrite = TRUE)
