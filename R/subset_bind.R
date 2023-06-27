#' subset_bind
#'
#' @inheritParams step_lead_lag
#' @inheritParams recipes::bake
#'
#' @param to_bind matrix to column bind
#'
#' @return columnwise combined tibble of `new_data` and `to_bind`
#'
#' @noRd
subset_bind <- function(new_data, to_bind, n_subset, n_shift) {
  if (n_subset > 1) {
    ind <- seq(
      n_shift + 1,
      nrow(new_data),
      n_subset
    )

    return(bind_cols(new_data[ind, ],
      as_tibble(to_bind, .name_repair = "minimal"),
      .name_repair = "minimal"
    ))
  }


  bind_cols(new_data, to_bind)
}
