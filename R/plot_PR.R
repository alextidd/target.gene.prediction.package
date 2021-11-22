plot_PR <- function(PR, ...){
  PR %>%
    dplyr::filter(.threshold %ni% c(0, Inf, -Inf),
                  recall != 1,
                  precision != 0) %>%
    dplyr::mutate(line = dplyr::n() > 1,
                  point = dplyr::n() == 1) %>%
    ggplot2::ggplot(ggplot2::aes(x = recall,
                                 y = precision,
                                 ...)) +
    ggplot2::geom_line(data = . %>% dplyr::filter(prediction_type == "score")) +
    ggplot2::geom_point(data = . %>% dplyr::filter(prediction_type %in% c("max", "prediction")),
                        ggplot2::aes(shape = prediction_type)) +
    ggplot2::xlim(0,1) +
    ggplot2::ylim(0,1) +
    ggplot2::coord_equal()
}


# PR %>%
#   dplyr::filter(.threshold %ni% c(0, Inf, -Inf),
#                 recall != 1,
#                 precision != 0) %>%
#   dplyr::mutate(line = dplyr::n() > 1,
#                 point = dplyr::n() == 1) %>%
#   ggplot2::ggplot(ggplot2::aes(x = recall,
#                                y = precision,
#                                ...)) +
#   ggplot2::geom_line(data = . %>% dplyr::filter(line)) +
#   ggplot2::geom_point(data = . %>% dplyr::filter(point)) +
#   ggplot2::xlim(0,1) +
#   ggplot2::ylim(0,1) +
#   ggplot2::coord_equal()