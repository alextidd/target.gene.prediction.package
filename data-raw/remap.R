## code to prepare `remap` dataset goes here

ReMapDir <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/ReMap/split_by_Experiment.TF.Biotype/"
remap <- list()
for(file in (list.files(ReMapDir, full.names = T) %>% head)) {
  name <- basename(file)
  remap[[name]] <- target.gene.prediction.package::import_BED_to_GRanges(file, metadata_cols = c("Experiment.TranscriptionFactor.CellType"))
  # Split metadata string into columns
  mcols(remap[[name]]) <- mcols(remap[[name]]) %>%
    dplyr::as_tibble() %>%
    tidyr::separate(col = Experiment.TranscriptionFactor.CellType,
                    into = c("Experiment", "TranscriptionFactor", "CellType"),
                    sep = "\\.")
}
# Unlist
remap <- unlist(as(remap, "GRangesList"))

usethis::use_data(remap, overwrite = TRUE)
