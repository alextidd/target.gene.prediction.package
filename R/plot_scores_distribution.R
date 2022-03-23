plot_scores_distribution <- function(df, score_col){
  df %>%
    tidyr::unite(pair, cs, enst) %>%
    dplyr::select(pair, score = {{score_col}}, known_gene) %>%
    dplyr::group_by(pair) %>%
    dplyr::filter(score == max(score)) %>%
    dplyr::ungroup() %>%
    dplyr::distinct() %>%
    ggplot2::ggplot(ggplot2::aes(x = reorder(pair, -score),
                                 y = score,
                                 fill = known_gene)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::geom_point(data = . %>% dplyr::filter(known_gene),
                        colour = "black") +
    ggplot2::labs(title = annotation) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
}
