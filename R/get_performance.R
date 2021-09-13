#' Get performance metrics for the predictions
#'
#' Calculate p-value, OR, precision, recall, sensitivity, specificity and F-score for the prediction.
#'
#' @param predictions Scores df in predict_target_genes function
#'
#' @return `performance` object (performance)
#' @export
get_performance <- function(predictions){
  data.frame(Positive = predictions %>% dplyr::filter(prediction) %>% dplyr::pull(symbol) %>% dplyr::n_distinct(),
             True = predictions %>% condition_n_genes(driver),
             TP = predictions %>% condition_n_genes(prediction & driver),
             FP = predictions %>% condition_n_genes(prediction & !driver),
             FN = predictions %>% condition_n_genes(!prediction & driver),
             TN = predictions %>% condition_n_genes(!prediction & !driver)) %>%
    dplyr::mutate(p = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$p.value,
                  OR = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$estimate,
                  Precision = TP / (TP + FN),
                  Recall = TP / (TP + FP),
                  Sensitivity = TP / (TP + FP),
                  Specificity = TN / (TN + FP),
                  Fscore = (Precision * Recall) / (Precision + Recall))
}

condition_n_genes <- function(df, ...){
  df %>% dplyr::filter(...) %>% dplyr::pull(symbol) %>% dplyr::n_distinct()
}
