## code to prepare `sysdata` dataset goes here
library(devtools) ; setwd("/working/lab_jonathb/alexandT/tgp") ; load_all()
reference_panels_dir <- "/working/lab_jonathb/alexandT/tgp_paper/wrangle_package_data/reference_panels/"
GENCODEPath <- paste0(reference_panels_dir, "/output/GENCODE/")
REVELPath <- paste0(reference_panels_dir, "/output/REVEL/")

# annotations_metadata
# manually copied from
# https://docs.google.com/spreadsheets/d/1De71B8qUdNge9jrH65GvryX6eYYHgdHy707HtDxyU2k/edit#gid=1075783341
# > data/metadata.tsv
annotations_metadata <- read_tibble("data/metadata.tsv", header = T)

# GENCODE
pcENSGs <- read_tibble(paste0(GENCODEPath, "proteincoding.gencode.v34lift37.basic.ENSGs.txt"))$V1
TSSs <- import_BED(gzfile(paste0(GENCODEPath, "gencode.v34lift37.basic.tss.bed.gz")), metadata_cols = c("ensg", "symbol", "enst"))
promoters <- import_BED(gzfile(paste0(GENCODEPath, "gencode.v34lift37.basic.promoter.bed.gz")), metadata_cols = "enst")
introns <- import_BED(gzfile(paste0(GENCODEPath, "gencode.v34lift37.basic.intron.bed.gz")), metadata_cols = "enst")
exons <- import_BED(gzfile(paste0(GENCODEPath, "gencode.v34lift37.basic.exon.bed.gz")), metadata_cols = "enst")

# ChrSizes
ChrSizes <- read_tibble(paste0(reference_panels_dir, "data/hg/hg19.genome"))
names(ChrSizes) <- c("chrom", "size")

# REVEL - deleterious coding variants
missense <- paste0(REVELPath, "missense.tsv") %>% read_tibble(header = T)
nonsense <- paste0(REVELPath, "nonsense.tsv") %>% read_tibble(header = T)
splicesite <- paste0(REVELPath, "splicesite.tsv") %>% read_tibble(header = T)

usethis::use_data(annotations_metadata,
                  pcENSGs,
                  TSSs,
                  promoters,
                  introns,
                  exons,
                  ChrSizes,
                  missense,
                  nonsense,
                  splicesite,
                  internal = TRUE,
                  overwrite = TRUE)
