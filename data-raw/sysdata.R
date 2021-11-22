## code to prepare `sysdata` dataset goes here
library(devtools) ; setwd("/working/lab_georgiat/alexandT/tgp") ; load_all()

# introns
intronFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/reference_data/GENCODE/proteincoding.gencode.v34lift37.basic.intron.bed.gz"
introns <- import_BED(gzfile(intronFile), metadata_cols = c("enst"))

# exons
exonFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/reference_data/GENCODE/proteincoding.gencode.v34lift37.basic.exon.bed.gz"
exons <- import_BED(gzfile(exonFile), metadata_cols = c("enst"))

# promoters
promFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/reference_data/GENCODE/proteincoding.gencode.v34lift37.basic.promoter.bed.gz"
promoters <- import_BED(gzfile(promFile),  metadata_cols = c("enst"))

# TSSs
TSSFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/reference_data/GENCODE/proteincoding.TSSs.GENCODE.bed.gz"
TSSs <- import_BED(gzfile(TSSFile), metadata_cols = c("ensg", "symbol", "enst"))

# ChrSizes
ChrSizesFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/reference_data/hg19.genome"
ChrSizes <- read_tibble(ChrSizesFile)
names(ChrSizes) <- c("chrom", "size")

usethis::use_data(introns,
                  exons,
                  promoters,
                  TSSs,
                  ChrSizes,
                  internal = TRUE)
