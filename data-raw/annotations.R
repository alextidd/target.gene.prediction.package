## code to prepare `annotations` dataset goes here

# A list of various annotations to be intersected with SNPs and genes for pair predictions

# DHS mark binning dataset
DHSsDir <- "/working/lab_georgiat/jonathB/PROJECTS/trench_lab/target_gene_prediction/output/activity_signal_matrix/bin_regions/"
# load
DHSs <- list() ; for(Method in c("specificity")) {
  for(Mark in c("H3K27ac")){
    dir <- paste0(DHSsDir, "/bin_regions_by_", Method, "/sorted/", Mark, "/quartiles/")
    for(file in list.files(dir, full.names = T)){
      name <- file %>% basename %>% sub(".sorted.bed","",.)
      DHSs[[name]] <- target.gene.prediction.package::import_BED(file)
      # Add metadata columns
      DHSs[[name]]$Method <- Method
      DHSs[[name]]$Mark <- Mark
      DHSs[[name]]$CellType <- sub("\\..*", "", name)
      DHSs[[name]]$Binning <- "quartiles"
      DHSs[[name]]$Bin <- sub(".*\\.bin", "", name) %>% as.numeric
      # TODO: fix data structure so that the metadata is within the bed files, not extracted from the path info
      # TODO: fix bin columns within the files (character -> numeric)
      # TODO: fix cell type naming codes
    }
  }
} ; DHSs <- dplyr::bind_rows(DHSs)

# TFBSs ReMap dataset (breast)
ReMapDir <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/ReMap/split_by_Experiment.TF.Biotype/"
ReMapFiles <- list.files(ReMapDir, full.names = T)
breast_ReMapFiles <- ReMapFiles[greplany(target.gene.prediction.package::annotations_metadata[["TFBSs"]] %>%
                                           dplyr::filter(tissue=="breast") %>%
                                           dplyr::pull(code),
                                         ReMapFiles)]
# load
TFBSs <- list() ; for(file in breast_ReMapFiles) {
  name <- basename(file)
  TFBSs[[name]] <- target.gene.prediction.package::import_BED(file,
                                                              metadata_cols = c("Experiment.TranscriptionFactor.CellType"))
  # Split metadata string into columns
  TFBSs[[name]] <- TFBSs[[name]] %>%
    tidyr::separate(col = Experiment.TranscriptionFactor.CellType,
                    into = c("Experiment", "TranscriptionFactor", "CellType"),
                    sep = "\\.",
                    extra = "merge")
} ; TFBSs <- dplyr::bind_rows(TFBSs)

# List together annotations
annotations <- list("DHSs" = DHSs,
                    "TFBSs" = TFBSs)

# Save
usethis::use_data(annotations, overwrite = TRUE)

