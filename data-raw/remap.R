## code to prepare `remap` dataset goes here

ReMapDir <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/ReMap/split_by_Experiment.TF.Biotype/"
remap <- list()
for(file in (list.files(ReMapDir, full.names = T) %>% head)) {
  name <- basename(file)
  remap[[name]] <- target.gene.prediction.package::import_BED(file,
                                                              metadata_cols = c("Experiment.TranscriptionFactor.CellType"))
  # Split metadata string into columns
  remap[[name]] <- remap[[name]] %>%
    tidyr::separate(col = Experiment.TranscriptionFactor.CellType,
                    into = c("Experiment", "TranscriptionFactor", "CellType"),
                    sep = "\\.")
}
# Unlist
remap <- bind_rows(remap)

usethis::use_data(remap, overwrite = TRUE)
