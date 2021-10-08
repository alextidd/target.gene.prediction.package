#' Get PR statistics for the predictions
#'
#' Get PR and AUPRC for predictions
#'
#' @param predictions `predictions` df in predict_target_genes function
#' @param ... Columns to use as predictors
#'
#' @return `PR` tibble
#' @export
get_PR <- function(scores, txv_master, ...){

  performance <- list()

  # score, prediction and max
  pred <- scores %>%
    # get predictions only in CSs with a driver within max prediction distance for performance evaluation
    get_testable() %>%
    dplyr::select(cs, symbol, ..., driver) %>%
    # gather prediction method columns
    tidyr::pivot_longer(names_to = "prediction_method",
                        values_to = "score",
                        ...) %>%
    # replace NAs with 0s
    dplyr::mutate(score = score %>% tidyr::replace_na(0)) %>%
    # get maximum score per CS-x-gene-x-method
    dplyr::group_by(prediction_method, cs, symbol) %>%
    dplyr::filter(score == max(score)) %>%
    # find maximum scoring gene per CS-x-method
    dplyr::group_by(prediction_method, cs) %>%
    dplyr::mutate(max = as.numeric(score == max(score))) %>%
    dplyr::ungroup() %>%
    # find positives (score > median(score)) overall
    dplyr::mutate(prediction = as.numeric(score > median(score))) %>%
    # gather prediction type columns
    tidyr::pivot_longer(c(score, max, prediction),
                        names_to = "prediction_type",
                        values_to = "prediction")

  # get summary statistics (various performance metrics)
  performance$summary <- pred %>%
    # score is not binary, cannot be summarised, filter out
    dplyr::filter(prediction_type != "score") %>%
    dplyr::mutate(prediction = as.logical(prediction)) %>%
    dplyr::group_by(prediction_method, prediction_type) %>%
    dplyr::group_modify(
      ~ data.frame(True = .x %>% condition_n_genes(driver),
                   Positive = .x %>% condition_n_genes(prediction),
                   TP = .x %>% condition_n_genes(prediction & driver),
                   FP = .x %>% condition_n_genes(prediction & !driver),
                   TN = .x %>% condition_n_genes(!prediction & !driver),
                   FN = .x %>% condition_n_genes(!prediction & driver))
    ) %>%
    dplyr::mutate(p = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$p.value,
                  OR = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$estimate,
                  Precision = TP / (TP + FN),
                  Recall = TP / (TP + FP),
                  Sensitivity = TP / (TP + FP),
                  Specificity = TN / (TN + FP),
                  Fscore = (Precision * Recall) / (Precision + Recall))

  # format for PR function input
  PR_in <- pred %>%
    dplyr::select(prediction_type, prediction_method, prediction, driver) %>%
    dplyr::group_by(prediction_method, prediction_type) %>%
    # refactor driver predictions for PR function
    dplyr::mutate(driver = ifelse(driver, "positive", "negative") %>% factor(c("positive", "negative")))

  performance$PR <- PR_in %>%
    # calculate PR curve
    yardstick::pr_curve(driver, prediction) %>%
    dplyr::left_join(PR_in %>%
                       # calculate AUPRC
                       yardstick::pr_auc(driver, prediction) %>%
                       dplyr::select(prediction_method, prediction_type,
                                     AUC = .estimate))

  return(performance)

}
