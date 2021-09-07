## code to prepare `DHSs` dataset goes here

# DHS mark binning dataset - top decile for each + H3K27ac only (! Temporary - file too big otherwise)
DHSsDir <- "/working/lab_georgiat/jonathB/PROJECTS/trench_lab/target_gene_prediction/output/activity_signal_matrix/bin_regions/"
# load
DHSs <- list() ; for(Method in c("specificity", "signal")) {
  for(Mark in c("H3K27ac")){
    dir <- paste0(DHSsDir, "/bin_regions_by_", Method, "/sorted/", Mark, "/deciles/")
    for(file in list.files(dir, pattern = "bin10", full.names = T)){
      df <- target.gene.prediction.package::import_BED(file)
      # Add metadata columns
      CellType <- sub("\\..*", "", basename(file))
      Bin <- basename(file) %>% sub(".sorted.bed", "", .) %>% sub(".*\\.", "", .)
      # TODO: fix data structure so that the metadata is within the bed files, not extracted from the path info
      # TODO: fix bin columns within the files (character -> numeric)
      # TODO: fix cell type naming codes
      DHSs[[Method]][[Mark]][[CellType]][[Bin]] <- df
    }
  }
}

# Save
usethis::use_data(DHSs, overwrite = TRUE)







# # TFBSs ReMap dataset (breast)
# ReMapDir <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/ReMap/split_by_Experiment.TF.Biotype/"
# ReMapFiles <- list.files(ReMapDir, full.names = T)
# breast_ReMapFiles <- ReMapFiles[greplany(target.gene.prediction.package::annotations_metadata[["TFBSs"]] %>%
#                                            dplyr::filter(tissue=="breast") %>%
#                                            dplyr::pull(code),
#                                          ReMapFiles)]
# # load
# TFBSs <- list() ; for(file in breast_ReMapFiles) {
#   name <- basename(file)
#   TFBSs[[name]] <- target.gene.prediction.package::import_BED(file,
#                                                               metadata_cols = c("Experiment.TranscriptionFactor.CellType"))
#   # Split metadata string into columns
#   TFBSs[[name]] <- TFBSs[[name]] %>%
#     tidyr::separate(col = Experiment.TranscriptionFactor.CellType,
#                     into = c("Experiment", "TranscriptionFactor", "CellType"),
#                     sep = "\\.",
#                     extra = "merge")
# } ; TFBSs <- dplyr::bind_rows(TFBSs)
