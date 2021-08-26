## code to prepare `drivers` dataset goes here
library(dplyr)

DriversDir <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/Traits/BC/KnownGenes/"
drivers <- list()
for(file in list.files(DriversDir, full.names = T, pattern = "drivers")) {
  name <- basename(file)
  drivers[[name]] <- target.gene.prediction.package::read_tibble(file) %>%
    select("symbol" = V1)
}
drivers <- bind_rows(drivers, .id = "DriversList")

usethis::use_data(drivers, overwrite = TRUE)
