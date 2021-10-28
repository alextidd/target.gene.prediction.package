# This function takes vectors of MAE assay names to add / multiply together and generates predictions from
# those components with different combinatins of weightings. It then measures the performance of each model
# and returns a tibble with the models' weights and the models' performance at predicting drivers as the target
# genes of variants in nearby CSs
# model building: for(W in model_weights){
#                   ( to_add[[1]] * W[[1]] + to_add[[2]] * W[[2]] + ... ) *
#                   ( to_multiply[[1]] * to_multiply[[2]] * ... )
#                 }
weight_and_score_manually <- function(MA,
                                      celltype_of_interest,
                                      txv_master,
                                      drivers,
                                      to_add = NULL,
                                      to_multiply = NULL,
                                      n_unique_manual_weights = NULL){

  # celltype_of_interest <- "BRST.MCF7.CNCR" ; n_unique_manual_weights = NULL

  # subset MA
  sub_MA <- MultiAssayExperiment::subsetByColData(MA, c("value", celltype_of_interest))
  pair_info <- txv_master %>% dplyr::select(pair, variant, cs, symbol)

  # get all possible weights
  model_weights <- to_add %>%
    sapply(function(x){
      seq(1, 100,
          length = ifelse(is.null(n_unique_manual_weights), length(to_add), n_unique_manual_weights))
    },
    simplify = F, USE.NAMES = T) %>%
    expand.grid %>%
    tibble::as_tibble()

  # score on each possible weight combination
  model_performance <- 1:nrow(model_weights) %>% #sample(10) %>%
    lapply(function(i) {
      print(i)
      # weight and score
      curr_weights <- model_weights %>% dplyr::filter(dplyr::row_number() == i)
      (
        if(is.null(to_add)){ 1 } else {
        to_add %>%
          lapply(function(a){ MultiAssayExperiment::assay(sub_MA, a) * curr_weights[[a]] }) %>%
          Reduce(`+`, .)
        }
      ) * (
        if(is.null(to_multiply)){ 1 } else {
          to_multiply %>%
            lapply(function(a){ MultiAssayExperiment::assay(sub_MA, a) }) %>%
            Reduce(`*`, .)
        }
      ) -> totals
      colnames(totals) <- "score"
      scores <- totals %>%
        tibble::as_tibble(rownames = "pair") %>%
        dplyr::left_join(pair_info, by = "pair") %>%
        dplyr::group_by(cs, symbol) %>%
        dplyr::summarise(score = max(score), .groups = "drop")

      # test performance
      out <- scores %>%
        get_PR(txv_master, drivers, score) %>%
        purrr::map(~ dplyr::mutate(., model = i))

      # return
      return(out)
    })

  # create list
  manual_models <- model_performance %>%
    purrr::map(~ dplyr::full_join(.$PR %>%
                                    dplyr::filter(prediction_type == "max") %>%
                                    dplyr::distinct(model, PR_AUC),
                                  .$summary %>%
                                    dplyr::filter(prediction_type == "max"),
                                  by = "model")) %>%
    dplyr::bind_rows() %>%
    dplyr::select(model, everything(), -prediction_method, -prediction_type) %>%
    dplyr::arrange(-PR_AUC) %>%
    dplyr::left_join(model_weights %>% dplyr::mutate(model = dplyr::row_number()), by = "model")

  return(manual_models)
}
