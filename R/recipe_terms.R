#' recipe_term_names
#'
#' @description Get the names of columns that are added grouped by the step.
#'
#' @inheritParams predict_terms
#' @param object either a 'recipe' or 'step' with grouped regressors.
#'
#' @details In some cases regressors should be treated like a group.  This
#' often occurs when dealing with lagged relationships, and when doing signal
#' decomposition. This function attempts to get the names of new columns
#' created as a part of each step and returns it as a list having the same
#' length as the number of steps.
#'
#' @return A list of character vectors grouped by step.
#'
#' @export
#' @examples
#' data(transducer)
#' transducer$datetime_num <- as.numeric(transducer$datetime)
#'
#' rec_toll_rasmussen <- recipe(wl~baro + datetime_num, transducer) |>
#'    step_lead_lag(baro, lag = log_lags(100, 86400 * 2 / 120)) |>
#'    prep()
#'
#' recipe_term_names(rec_toll_rasmussen)
recipe_term_names <- function(object, ...) UseMethod("recipe_term_names")

#' @rdname recipe_term_names
#' @export
recipe_term_names.recipe <- function(object, ...) {
  lapply(object$steps, function(x) {
    recipe_term_names(x)
  })
}

#' @rdname recipe_term_names
#' @export
recipe_term_names.step_lag <- function(object, ...) {

  columns <- object$columns
  nms     <- vector(mode = "list", length = length(columns))

  for(i in seq_along(columns)) {
    nms[[i]] <- paste0(object$prefix, object$lag, '_', columns[i])
  }

  return(nms)
}

#' @rdname recipe_term_names
#' @export
recipe_term_names.step_distributed_lag <- function(object, ...) {
  new_names <- list()

  for (i in seq_along(object$columns)) {

    new_names[[i]] <- paste0(object$prefix, object$columns[i], '_',
                             1:ncol(object$basis_mat))
  }

  return(new_names)

}

#' @rdname recipe_term_names
#' @export
recipe_term_names.step_lead_lag <- function(object, ...) {

  return(recipe_term_names.step_lag(object))

}


#' @rdname recipe_term_names
#' @export
recipe_term_names.step_harmonic <- function(object, ...) {
  col_names <- names(object$starting_val)
  n_frequency <- length(object$frequency)

  for (i in seq_along(col_names)) {

    new_names[[i]] <- paste0(
      col_names[i],
      rep(c("_sin_", "_cos_"), each = n_frequency),
      1:n_frequency
    )
  }
  return(new_names)

}

#' @rdname recipe_term_names
#' @export
recipe_term_names.step_earthtide <- function(object, ...) {
  n_terms <- nrow(object$wave_groups)

  new_names <- paste0(object$prefix,
                      rep(c('sin_', 'cos_'), each = n_terms),
                      rep(1:n_terms, times = 2))
  return(list(new_names))

}

#' @rdname recipe_term_names
#' @export
recipe_term_names.step_ns <- function(object, ...) {

  new_names <- list()
  for (i in 1:length(object$objects)) {
    orig_var <- attr(object$objects[[i]], "var")
    new_cols <- vapply(object$objects, ncol, c(int = 1L))
    new_names[[i]] <-
      paste(orig_var, "ns", names0(new_cols[i], ""), sep = "_")
  }

  return(new_names)

}



