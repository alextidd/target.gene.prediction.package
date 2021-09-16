plot_PR <- function(PR){
  PR %>%
    dplyr::filter(.threshold %ni% c(0, Inf, -Inf) & recall != 1) %>%
    dplyr::mutate(line = dplyr::n() > 1,
                  point = dplyr::n() == 1) %>%
    ggplot2::ggplot(ggplot2::aes(x = precision,
                                 y = recall,
                                 colour = prediction_type)) +
    ggplot2::geom_line(data = . %>% dplyr::filter(line)) +
    ggplot2::geom_point(data = . %>% dplyr::filter(point))
}
