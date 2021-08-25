#' Import a BEDPE
#'
#' Import a paired-end BEDPE file as a list object c(dataframe of first end, dataframe of last end)
#'
#' @param bedfile A tab-delimited file in bedpe format with columns: chrA, startA, stopA, chrB, startB, stopB, ...(metadata cols)...
#' @param metadata_cols A character vector of the names of metadata columns (bedfile cols 4+) in the BED file
#'
#' @return A paired (first + last) list object representation of the BEDPE file with an ID column matching pairs, optionally with metadata columns
#' @export
#' @importFrom dplyr %>%
#' @importFrom rlang .data
import_BEDPE_to_List <- function(bedpefile, metadata_cols = NULL) {
  message(paste("Importing paired-end BED file from", bedpefile))

  bedpe <- target.gene.prediction.package::read_tibble(bedpefile) %>%
    mutate(PairID = dplyr::row_number())
  PEList <- list()
  PEList[["first"]] <- bedpe[,-c(4:6)]
  PEList[["last"]] <- bedpe[,-c(1:3)]
  lapply(PEList,
         function(x)  x %>%
           dplyr::rename_with(.data, ~ c("chr", "start", "end"), 1:3) %>%
           { if(!is.null(metadata_cols)) dplyr::rename_with(.data, ~ metadata_cols, 4:(3 + length(metadata_cols))) else .data }
  )
}
