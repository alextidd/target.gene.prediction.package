## code to prepare `introns` dataset goes here

intronFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/reference_data/GENCODE/proteincoding.gencode.v34lift37.basic.intron.bed.gz"
introns <- target.gene.prediction.package::import_BED(gzfile(intronFile),
                                                        metadata_cols = c("enst"))

usethis::use_data(introns, overwrite = TRUE)
