## code to prepare `sysdata` dataset goes here
library(devtools) ; setwd("/working/lab_georgiat/alexandT/tgp") ; load_all()
GENCODEPath="/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/reference_data/GENCODE/proteincoding.gencode.v34lift37.basic."

# fix up GENCODE annotations files for compatibility
wrangle_GENCODE <- function(ann, keepAllMetadata = F){
  ann %>%
    # convert 1-based -> 0-based
    dplyr::mutate(start = start - 1) %>%
    # remove version number extensions from ensts
    { if("enst" %in% colnames(ann)) dplyr::mutate(., enst = enst %>% tools::file_path_sans_ext()) else . } %>%
    # remove ENSG and symbol columns unless specified
    { if(!keepAllMetadata) dplyr::select(., -c("ensg", "symbol")) else . }
}

# TSSs
TSSFile <- paste0(GENCODEPath, "tss.bed.gz")
TSSs <- import_BED(gzfile(TSSFile), metadata_cols = c("ensg", "symbol", "enst")) %>%
  wrangle_GENCODE(keepAllMetadata = T)

# promoters
promFile <- paste0(GENCODEPath, "promoter.bed.gz")
promoters <- import_BED(gzfile(promFile),  metadata_cols = c("ensg", "symbol", "enst")) %>%
  wrangle_GENCODE

# introns
intronFile <- paste0(GENCODEPath, "intron.bed.gz")
introns <- import_BED(gzfile(intronFile), metadata_cols = c("ensg", "symbol", "enst")) %>%
  wrangle_GENCODE

# exons
exonFile <- paste0(GENCODEPath, "exon.bed.gz")
exons <- import_BED(gzfile(exonFile), metadata_cols = c("ensg", "symbol", "enst")) %>%
  wrangle_GENCODE

# ChrSizes
ChrSizesFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/reference_data/hg19.genome"
ChrSizes <- read_tibble(ChrSizesFile)
names(ChrSizes) <- c("chrom", "size")

# REVEL - deleterious coding variants - only those found in GENCODE Basic protein-coding # TODO: fix this in upstream scripts
REVELPath="/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/REVEL/"
missenseFile <- paste0(REVELPath, "missense_SNVs_dbNSFPv4.1a.tsv")
nonsenseFile <- paste0(REVELPath, "nonsense_SNVs_dbNSFPv4.1a.tsv")
splicesiteFile <- paste0(REVELPath, "splicesite_SNVs.tsv")

wrangle_REVEL <- function(df){
  df %>%
    tidyr::separate_rows(ensgs) %>%
    dplyr::group_by(dplyr::across(c(-ensgs))) %>%
    dplyr::summarise(ensgs = ensgs %>% unique %>% paste(collapse = ";")) %>%
    dplyr::ungroup()
}

missense <- missenseFile %>% read_tibble %>%
  dplyr::select(chrom = V1, position = V2, ensgs = V3, score = V4) %>%
  wrangle_REVEL
nonsense <- nonsenseFile %>% read_tibble %>%
  dplyr::select(chrom = V1, position = V2, ensgs = V3) %>%
  wrangle_REVEL
splicesite <- splicesiteFile %>% read_tibble %>%
  dplyr::select(chrom = V1, position = V2, ensgs = V3) %>%
  wrangle_REVEL

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
