## code to prepare `promoters` dataset goes here

promFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/reference_data/GENCODE/proteincoding.gencode.v34lift37.basic.promoter.bed.gz"
promoters <- target.gene.prediction.package::import_BED(gzfile(promFile),
                                                        metadata_cols = c("enst"))

usethis::use_data(promoters, overwrite = TRUE)
