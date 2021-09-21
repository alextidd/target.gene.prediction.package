#' Get performance metrics for the predictions
#'
#' Limit predictions to credible sets with at least one driver within the maximum prediction distance.
#' Calculate p-value, OR, precision, recall, sensitivity, specificity and F-score for the prediction.
#'
#' @param predictions Scores df in predict_target_genes function
#'
#' @return `performance` object (performance)
#' @export
get_performance <- function(predictions, prediction_col){

  # Simplify to CS-gene interactions and get maximum score for each
  predictions <- predictions %>%
    dplyr::group_by(cs, symbol) %>%
    dplyr::filter(score == max(score))

  # Calculate p-value, OR, precision, recall, sensitivity, specificity and F-score for the prediction
  data.frame(
    Positive = predictions %>% condition_n_genes({{prediction_col}}),
    True = predictions %>% condition_n_genes(driver),
    TP = predictions %>% condition_n_genes({{prediction_col}} & driver),
    FP = predictions %>% condition_n_genes({{prediction_col}} & !driver),
    TN = predictions %>% condition_n_genes(!{{prediction_col}} & !driver),
    FN = predictions %>% condition_n_genes(!{{prediction_col}} & driver)) %>%
    dplyr::mutate(p = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$p.value,
                  OR = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$estimate,
                  Precision = TP / (TP + FN),
                  Recall = TP / (TP + FP),
                  Sensitivity = TP / (TP + FP),
                  Specificity = TN / (TN + FP),
                  Fscore = (Precision * Recall) / (Precision + Recall))

}

