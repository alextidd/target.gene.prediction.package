get_testable <- function(df){
  df %>%
    dplyr::group_by(cs) %>%
    dplyr::filter(any(driver==T)) %>%
    dplyr::ungroup()
}
