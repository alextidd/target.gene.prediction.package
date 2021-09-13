#' Get PR statistics for the predictions
#'
#' Get PR and AUPRC for predictions
#'
#' @param predictions `predictions` df in predict_target_genes function
#' @param ... Columns to use as predictors
#'
#' @return `PR` tibble
#' @export
get_PR <- function(predictions, ...){

  PR_in <- predictions %>%
    dplyr::mutate(driver = ifelse(driver, "positive", "negative") %>% factor(c("positive", "negative"))) %>%
    # if score is logical, convert to binary
    dplyr::mutate_if(is.logical, as.numeric) %>%
    # gather and group
    dplyr::select(..., driver) %>%
    tidyr::gather(key = "prediction_type",
                  value = "prediction",
                  ...) %>%
    dplyr::group_by(prediction_type)

  PR <- PR_in %>%
    # calculate PR curve
    yardstick::pr_curve(driver, prediction) %>%
    dplyr::left_join(PR_in %>%
                       # calculate AUPRC
                       yardstick::pr_auc(driver, prediction) %>%
                       dplyr::select(prediction_type,
                                     AUC = .estimate))

  return(PR)

}
