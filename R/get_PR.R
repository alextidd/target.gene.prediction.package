#' Get PR statistics for the predictions
#'
#' Get PR and AUPRC for predictions
#'
#' @param scores scores object passed from predict_target_genes
#' @param txv_master txv_master object passed from predict_target_genes
#' @param drivers drivers object passed from predict_target_genes
#'
#' @return `PR` tibble
#' @export
get_PR <- function(scores, txv_master, drivers){

  performance <- list()

  # get all testable CS-gene pairs
  testable <- txv_master %>%
    # add drivers
    dplyr::mutate(driver = symbol %in% drivers$symbol) %>%
    # get predictions only in CSs with a driver within max prediction distance for performance evaluation
    get_testable() %>%
    dplyr::distinct(cs, symbol, driver)

  # score, max, max_score
  pred <- scores %>%
    dplyr::select(-c(chrom:end)) %>%
    dplyr::right_join(testable) %>%
    dplyr::group_by(cs, symbol, driver) %>%
    # get maximum score per CS-x-gene-x-method
    dplyr::summarise(
      dplyr::across(where(is.numeric), ~ max(.x))
    ) %>%
    # gather prediction methods
    tidyr::pivot_longer(
      where(is.numeric),
      names_to = "prediction_method",
      values_to = "score"
    ) %>%
    # max prediction
    dplyr::group_by(prediction_method, cs) %>%
    dplyr::mutate(max = as.numeric(score == max(score) & score > 0),
                  max_score = dplyr::case_when(max == 0 ~ 0,
                                               TRUE ~ score)) %>%
    # gather prediction types
    tidyr::pivot_longer(
      c(score, max_score, max),
      names_to = "prediction_type",
      values_to = "prediction"
    ) %>%
    dplyr::ungroup()

  # get summary statistics (various performance metrics)
  performance$summary <- pred %>%
    # score is not binary, cannot be summarised, filter out
    dplyr::filter(!grepl("score", prediction_type)) %>%
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

  # Add area under curve metric to summary
  performance$summary <- performance$summary %>%
    dplyr::left_join(performance$PR %>% dplyr::distinct(prediction_method,
                                                        prediction_type,
                                                        PR_AUC))
  return(performance)

}





