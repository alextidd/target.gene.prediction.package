## code to prepare `sysdata` dataset goes here
library(devtools) ; setwd("/working/lab_jonathb/alexandT/tgp") ; load_all()
GENCODEPath="/working/lab_jonathb/alexandT/tgp_paper/output/GENCODE/proteincoding.gencode.v34lift37.basic."
REVELPath="/working/lab_jonathb/alexandT/tgp_paper/output/REVEL/"

# GENCODE
TSSs <- import_BED(gzfile(paste0(GENCODEPath, "tss.bed.gz")), metadata_cols = c("ensg", "symbol", "enst"))
promoters <- import_BED(gzfile(paste0(GENCODEPath, "promoter.bed.gz") ), metadata_cols = "enst")
introns <- import_BED(gzfile(paste0(GENCODEPath, "intron.bed.gz")), metadata_cols = "enst")
exons <- import_BED(gzfile(paste0(GENCODEPath, "exon.bed.gz")), metadata_cols = "enst")

# ChrSizes
ChrSizes <- read_tibble( "/working/lab_jonathb/alexandT/tgp_paper/data/hg19.genome")
names(ChrSizes) <- c("chrom", "size")

# REVEL - deleterious coding variants - only those found in GENCODE Basic protein-coding # TODO: fix this in upstream scripts
missense <- paste0(REVELPath, "missense_SNVs.tsv") %>% read_tibble(header = T)
nonsense <- paste0(REVELPath, "nonsense_SNVs.tsv") %>% read_tibble(header = T)
splicesite <- paste0(REVELPath, "splicesite_SNVs.tsv") %>% read_tibble(header = T)

usethis::use_data(introns,
                  exons,
                  promoters,
                  TSSs,
                  ChrSizes,
                  missense,
                  nonsense,
                  splicesite,
                  internal = TRUE,
                  overwrite = TRUE)
