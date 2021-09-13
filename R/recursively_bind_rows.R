#' Recursively bind rows
#'
#' Recursively bind rows of a nested list.
#'
#' @param nested_list The nested list. The list must have a uniform nesting structure (every branch goes to the same depth).
#' @param nest_names The names of each successive nesting of the list.
#' @param unnest_dataframes Logical. Data frames are essentially lists. The unnesting will preserve data frames in the deepest nest unless this is set to TRUE.
#'
#' @return A recursively row-bound tibble, with nested list names as added ID columns.
#' @export
recursively_bind_rows <- function(nested_list, nest_names, unnest_dataframes = FALSE){

  # check that list is provided
  if(!is.list(nested_list)){ stop("Non-list supplied") }
  if(is.data.frame(nested_list) & unnest_dataframes == FALSE){ return(nested_list) }

  # internal function to find maximum depth of a nested list of dataframes
  depth_nested_list <- function(this, unnest_dfs = unnest_dataframes){
    if(is.data.frame(this) & !unnest_dfs){ # stop unnesting at second lowest level
      0L
    } else if(is.list(this)){
      1L + max(sapply(this, depth_nested_list))
    } else {
      0L
    }
  }

  # find max depth of the nesting and check that sufficient nest names are provided
  depth <- depth_nested_list(nested_list)
  if(length(nest_names) != depth){ stop("Number of nest names provided not equal to number of nests. ", depth, " names required.") }

  # start binding...
  unnested_list <- nested_list

  # bind deeper levels
  for(nest in (depth-1):1){
    # bind_rows applies to lists, so apply 1 level above the current deepest level
    nest_name = nest_names[[nest+1]]
    cat("depth =", nest+1, "\n")
    cat("nest name =", nest_name, "\n")
    unnested_list <- unnested_list %>%
      purrr::modify_depth(.depth = nest,
                          dplyr::bind_rows,
                          .id = nest_name)
  }

  # bind top nest
  cat("depth = 1\n")
  cat("nest name =",nest_names[1],"\n")
  unnested_list <- unnested_list %>%
    dplyr::bind_rows(.id = nest_names[1])

  # return
  return(unnested_list)

}
