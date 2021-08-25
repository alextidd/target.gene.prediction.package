#' Import a BED
#'
#' Import a BED file as a GRanges object
#'
#' @param bedfile A tab-delimited file in bed format with columns: chr, start, stop, ...(metadata cols)...
#' @param metadata_cols A character vector of the names of metadata columns (bedfile cols 4+) in the BED file, to be imported as metadata in the GRanges object
#'
#' @return A GRanges object representation of the BED file, optionally with metadata columns
#' @export
#' @importFrom rlang .data
import_BED_to_GRanges <- function(bedfile, metadata_cols = NULL) {
  message(paste("Importing BED file from", bedfile))

  target.gene.prediction.package::read_tibble(bedfile) %>%
    dplyr::rename(seqnames = .data$V1,
                  start = .data$V2,
                  end = .data$V3) %>%
    { if(!is.null(metadata_cols)) dplyr::rename_with(.data, ~ metadata_cols, 4:(3 + length(metadata_cols))) else .data } %>%
    plyranges::as_granges()
}
