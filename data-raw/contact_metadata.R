## code to prepare `contact_metadata` dataset goes here

library(dplyr)
contactDir="/working/lab_georgiat/alexandT/target_gene_prediction_paper/output/3C/"

contact_metadata <- dplyr::tibble(
  file = list.files(contactDir)
) %>%
  tidyr::separate(file,
                  sep = "\\_",
                  into = c("Study", "CellType", "Assay"),
                  remove = F) %>%
  mutate(Assay = Assay %>% sub("\\.bedpe$", "", .),
         contact_list_element = paste(Study, CellType, Assay, sep = "_")) %>%
# tissues
  mutate(Tissue = case_when(Study == "Javierre2016" ~ "blood",
                            Study %in% c("Orlando2018", "Chen2021") ~ "colorectal",
                            (Study == "Shi2021" & grepl("HaCaT", CellType)) ~ "epidermis",
                            (Study == "Shi2021" & CellType == "MyLa") ~ "blood",
                            CellType == "MCF7" ~ "breast"))

use_data(contact_metadata, overwrite = T)
