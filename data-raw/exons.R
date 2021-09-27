## code to prepare `exons` dataset goes here

exonFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/reference_data/GENCODE/proteincoding.gencode.v34lift37.basic.exon.bed.gz"
exons <- target.gene.prediction.package::import_BED(gzfile(exonFile),
                                                      metadata_cols = c("enst"))


usethis::use_data(exons, overwrite = TRUE)
