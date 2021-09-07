#' Recursively bind rows
#'
#' Recursively bind rows of a nested list.
#'
#' @param nested_list The nested list
#' @param nest_names The names of each successive nesting of the list. If none provided, nest1, nest2, nest3, ... will be used
#'
#' @return A recursively row-bound tibble, with nested list names as added ID columns.
#' @export
recursively_bind_rows <- function(nested_list, nest_names=NULL){

  # check that list is provided
  if(!is.list(nested_list)){ stop("Non-list supplied") }

  # internal function to find maximum depth of a nested list of dataframes
  depth_nested_dfs <- function(this){
    depth <- function(this) ifelse(is.list(this), 1L + max(sapply(this, depth)), 0L)
    depth(this) - 2
  }

  # find max depth of the nesting and check that sufficient nest names are provided
  depth <- depth_nested_dfs(nested_list)
  if(is.null(nest_names)){ nest_names = paste0("level", 1:(depth+1)) }
  if(length(nest_names) != depth+1){ stop("Number of nest names provided not equal to number of nests. ", depth+1, " names required.") }

  # start binding...
  unnested_list <- nested_list

  # bind deeper levels
  if(depth > 1){
    for(level in depth:1){
      nest_name = nest_names[[level+1]]
      cat("depth =", level+1, "\n")
      cat("nest name =", nest_name, "\n")
      unnested_list <- unnested_list %>%
        purrr::modify_depth(level, dplyr::bind_rows, .id = nest_name)
    }
  }

  # bind top nest
  cat("depth =", 1, "\n")
  cat("nest name =",nest_names[1],"\n")
  unnested_list <- unnested_list %>%
    dplyr::bind_rows(.id = nest_names[1])

  # return
  return(unnested_list)

}
