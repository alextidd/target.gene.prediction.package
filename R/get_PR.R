#' Get PR statistics for the predictions
#'
#' Get PR and AUPRC for predictions
#'
#' @param scores scores object passed from predict_target_genes
#' @param txv_master txv_master object passed from predict_target_genes
#' @param drivers drivers object passed from predict_target_genes
#' @param ... Columns to use as predictors (e.g. score column, values of individual annotation columns, ...)
#'
#' @return `PR` tibble
#' @export
get_PR <- function(scores, txv_master, drivers, ...){

  performance <- list()

  # score, prediction and max
  pred <- scores %>%
    # add drivers
    dplyr::mutate(driver = symbol %in% drivers$symbol) %>%
    # get predictions only in CSs with a driver within max prediction distance for performance evaluation
    get_testable() %>%
    dplyr::select(cs, symbol, ..., driver) %>%
    dplyr::distinct() %>%
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
    dplyr::mutate(max = as.numeric((score == max(score) & score != 0))) %>%
    dplyr::ungroup() %>%
    # find positives (score > median(score)) overall
    dplyr::mutate(prediction = as.numeric((score >= median(score) & score != 0))) %>%
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
      ~ data.frame(True = .x %>% condition_n_gene_x_cs_pairs(driver),
                   Positive = .x %>% condition_n_gene_x_cs_pairs(prediction),
                   TP = .x %>% condition_n_gene_x_cs_pairs(prediction & driver),
                   FP = .x %>% condition_n_gene_x_cs_pairs(prediction & !driver),
                   TN = .x %>% condition_n_gene_x_cs_pairs(!prediction & !driver),
                   FN = .x %>% condition_n_gene_x_cs_pairs(!prediction & driver))
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(p = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$p.value,
                  OR = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$estimate,
                  Precision = TP / (TP + FP),
                  Recall = TP / (TP + FN),
                  Sensitivity = TP / (TP + FN),
                  Specificity = TN / (TN + FP),
                  Fscore = (Precision * Recall) / (Precision + Recall)) %>%
    dplyr::ungroup()

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
                                     PR_AUC = .estimate)) %>%
    dplyr::ungroup()

  return(performance)

}
