## code to prepare `ChrSizes` dataset goes here

ChrSizesFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/reference_data/hg19.genome"
ChrSizes <- target.gene.prediction.package::read_tibble(ChrSizesFile)
names(ChrSizes) <- c("chrom", "size")

usethis::use_data(ChrSizes, overwrite = TRUE)
