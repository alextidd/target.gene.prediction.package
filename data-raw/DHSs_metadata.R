## code to prepare `annotations_metadata` dataset goes here

# Metadata for the 'annotations' list

# Metadata - H3K27ac-specificity-binned DHSs dataset
DHSs_metadata <- target.gene.prediction.package::read_tibble(
  "/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/DHS/metadata/DHS_corrected_biotypes.tsv",
  header = T) %>%
  dplyr::select(code,
                name,
                tissue,
                stage)


usethis::use_data(DHSs_metadata, overwrite = TRUE)


# # ReMap TFBSs dataset
# TFBSs <- target.gene.prediction.package::read_tibble(
#   "/working/lab_georgiat/alexandT/target_gene_prediction_paper/data/ReMap/metadata/remap_corrected_biotypes.tsv",
#   header = T) %>%
#   dplyr::select(code = biotype,
#                 name,
#                 tissue,
#                 stage)
