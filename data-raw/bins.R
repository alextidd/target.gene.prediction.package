## code to prepare `bins` dataset goes here

# DHS mark binning example dataset
BinsDir <- "/working/lab_georgiat/jonathB/PROJECTS/trench_lab/target_gene_prediction/output/activity_signal_matrix/bin_regions/bin_regions_by_specificity/sorted/H3K27ac/deciles/"
bins <- list()
for(file in list.files(BinsDir, pattern = "ENCFF", full.names = T)) {
  name <- file %>% basename %>% sub(".sorted.bed","",.)
  bins[[name]] <- target.gene.prediction.package::import_BED(file)
  # Split metadata string into columns
  bins[[name]]$Experiment.TranscriptionFactor.CellType <- name
  bins[[name]] <- bins[[name]] %>%
    tidyr::separate(col = Experiment.TranscriptionFactor.CellType,
                    into = c("CellType", "Mark", "Bin"),
                    sep = "\\.")
}
# Unlist
bins <- bind_rows(bins)

# (Temporary)
bins$Method <- "specificity"
bins$Binning <- "deciles"
bins$Bin <- bins$Bin %>% sub("bin","",.) %>% as.numeric
# TODO: fix data structure so that the metadata is within the bed files, not extracted from the path name
# TODO: fix bin columns within the files (character -> numeric)
# TODO: fix cell type naming codes

usethis::use_data(bins, overwrite = TRUE)
