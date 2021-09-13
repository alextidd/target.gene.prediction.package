## code to prepare `drivers` dataset goes here
library(dplyr)

DriversFile <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/Traits/BC/KnownGenes/breast_cancer_drivers_2021.txt"
drivers <- target.gene.prediction.package::read_tibble(DriversFile)$V1

usethis::use_data(drivers, overwrite = TRUE)


# DriversDir <- "/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/Traits/BC/KnownGenes/"
# drivers <- list()
# for(file in list.files(DriversDir, full.names = T, pattern = "drivers")) {
#   name <- basename(file)
#   drivers[[name]] <- target.gene.prediction.package::read_tibble(file) %>%
#     select("symbol" = V1)
# }
# drivers <- bind_rows(drivers, .id = "DriversList")
