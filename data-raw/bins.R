## code to prepare `bins` dataset goes here

# DHS mark binning example dataset
BinsDir <- "/working/lab_georgiat/jonathB/PROJECTS/trench_lab/target_gene_prediction/output/activity_signal_matrix/bin_regions/bin_regions_by_specificity/sorted/H3K27ac/deciles/"
bins <- list()
for(file in list.files(BinsDir, pattern = "ENCF", full.names = T)) {
  print(file)
  name <- file %>% basename %>% sub(".sorted.bed","",.)
  bins[[name]] <- target.gene.prediction.package::import_BED_to_GRanges(file)
  # Split metadata string into columns
  mcols(bins[[name]])$Experiment.TranscriptionFactor.CellType <- name
  mcols(bins[[name]]) <- mcols(bins[[name]]) %>%
    dplyr::as_tibble() %>%
    tidyr::separate(col = Experiment.TranscriptionFactor.CellType,
                    into = c("CellType", "Mark", "Bin"),
                    sep = "\\.")
}
# Unlist
bins <- unlist(as(bins, "GRangesList"))

# (Temporary)
bins$Method <- "specificity"
bins$Binning <- "deciles"
bins$Bin <- bins$Bin %>% sub("bin","",.) %>% as.numeric
# TODO: fix data structure so that the metadata is within the bed files, not extracted from the path name
# TODO: fix bin columns (character -> numeric)
# TODO: fix cell type naming codes

usethis::use_data(bins, overwrite = TRUE)
