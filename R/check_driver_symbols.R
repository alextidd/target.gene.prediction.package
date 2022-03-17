#' check that all symbols are in the GENCODE data
#'
#' @param drivers vector of driver symbols
#' @param driversFile name of the drivers file
#'
#' @return A file of variant-gene pair predictions, with associated scores, saved in the given output directory.
#' @export
check_driver_symbols <- function(drivers, driversFile){
  drivers_NotFound <- drivers[drivers %ni% TSSs$symbol]
  drivers_NonCoding <- dplyr::filter(TSSs,
                                     symbol %in% drivers[drivers %in% TSSs$symbol],
                                     ensg %ni% pcENSGs
                                     )$symbol
  if(length(c(drivers_NotFound, drivers_NonCoding)) > 0){
    message(dplyr::n_distinct(drivers_NotFound),
            " provided driver gene(s) from '", basename(driversFile), "' are not valid GENCODE gene symbols and therefore cannot be used.",
            "\nUnknown genes: ", paste(drivers_NotFound, sep="", collapse=", "))
    message(dplyr::n_distinct(drivers_NonCoding),
            " provided driver gene(s) from '", basename(driversFile), "' are GENCODE genes but are not protein-coding and therefore cannot be used.",
            "\nNon-coding genes: ", paste(drivers_NonCoding, sep="", collapse=", "))
    message(length(drivers[drivers %ni% c(drivers_NotFound, drivers_NonCoding)]),
            " drivers left after removing unknown/non-coding.")
  }

  TSSs %>%
    dplyr::filter(symbol %in% drivers, ensg %in% pcENSGs)
}
