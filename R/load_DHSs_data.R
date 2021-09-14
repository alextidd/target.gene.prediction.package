#' Load DHSs data
#'
#' Load the directory of DHSs data, downloaded from an FTP(?) alongside the package.
#' The function will take all files in this directory ending in *.bed, convert them to
#' BED tibble objects and add metadata.
#'
#' @param referenceDir The directory in which the package's accompanying data is stored.
#' (DHSs data is in the DHSs/ subdirectory).
#' Files must have columns ( chr | start | end | decile ) and *.bed file extension.
#'
#' @return `DHSs` object (list of BED tibbles containing binned DHSs)
#' @export
load_DHSs_data <- function(referenceDir){
  DHSsDir <- paste0(referenceDir, "/DHSs/")
  DHSs <- data.frame()
  for(file in list.files(DHSsDir, pattern = ".bed", full.names = T)){
    info <- stringr::str_split(basename(file), pattern = "\\.")[[1]]
    df <- target.gene.prediction.package::import_BED(file, metadata_cols = "Decile") %>%
      # Add metadata columns
      dplyr::mutate(Method = info[1],
                    Mark = info[2],
                    CellType = info[3])
    DHSs <- dplyr::bind_rows(DHSs, df)
  }
}
#
# DHSs = target.gene.prediction.package::DHSs %>%
#   target.gene.prediction.package::recursively_bind_rows(nest_names = c("Method", "Mark", "CellType", "Bin")) %>%
#   dplyr::mutate(annotation.description = paste(Bin, "of DHSs binned by", Method, "of", Mark, "annotations in",
#                                                target.gene.prediction.package::DHSs_metadata$name[target.gene.prediction.package::DHSs_metadata==CellType])
