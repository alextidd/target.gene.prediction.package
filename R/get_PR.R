get_PR <- function(scores, vxt_master, drivers, pcENSGs, max_n_drivers_per_CS){

  performance <- list()

  # get all testable CS-gene pairs
  testable <- vxt_master %>%
    # only test protein-coding target predictions against drivers (assumes all drivers are protein-coding)
    dplyr::filter(ensg %in% pcENSGs) %>%
    # add drivers
    dplyr::mutate(driver = symbol %in% drivers$symbol) %>%
    # only test predictions in CSs with a driver within max prediction distance for performance evaluation
    get_testable(max_n_drivers_per_CS) %>%
    dplyr::distinct(cs, symbol, driver)

  # score, max
  pred <- scores %>%
    dplyr::select(-dplyr::any_of(c("chrom", "start", "end"))) %>%
    dplyr::right_join(testable, by = c("cs", "symbol")) %>%
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
    dplyr::mutate(max = as.numeric(score == max(score) & score > 0)
                  #, max_score = dplyr::case_when(max == 0 ~ 0,
                  #                              TRUE ~ score)
                  ) %>%
    # gather prediction types
    tidyr::pivot_longer(
      c(score, max), #, max_score
      names_to = "prediction_type",
      values_to = "prediction"
    ) %>%
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
                       dplyr::select(prediction_method,
                                     prediction_type,
                                     PR_AUC = .estimate),
                     by = c("prediction_method", "prediction_type")) %>%
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
                   FN = .x %>% condition_n_gene_x_cs_pairs(!prediction & driver),
                   n_drivers = dplyr::filter(.x, driver)$symbol %>% dplyr::n_distinct())
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(p_value = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$p.value,
                  OR = fisher.test(matrix(c(TP,FP,FN,TN),2,2),alternative="greater")$estimate,
                  Precision = TP / (TP + FP),
                  Recall = TP / (TP + FN),
                  Sensitivity = TP / (TP + FN),
                  Specificity = TN / (TN + FP),
                  Fscore = (Precision * Recall) / (Precision + Recall)) %>%
    dplyr::ungroup() %>%
    # add area under curve metric to summary
    dplyr::left_join(performance$PR %>%
                       dplyr::distinct(prediction_method,
                                       prediction_type,
                                       PR_AUC),
                     by = c("prediction_method", "prediction_type"))

  return(performance)

}





