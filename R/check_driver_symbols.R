#' check that all symbols are in the GENCODE data
#'
#' @param drivers vector of driver symbols
#' @param driversFile name of the drivers file
#'
#' @return A file of variant-gene pair predictions, with associated scores, saved in the given output directory.
#' @export
check_driver_symbols <- function(drivers, driversFile){
  drivers_NotFound <- drivers[drivers %ni% target.gene.prediction.package::TSSs$symbol]
  if(length(drivers_NotFound) > 0){
    message(dplyr::n_distinct(drivers_NotFound),
            " provided driver gene(s) from '", basename(driversFile), "' are not found in the GENCODE gene symbols and therefore cannot be considered.",
            "\nUnknown genes: ", paste(drivers_NotFound, sep="", collapse=", "))
  }
}
