#' Create a lead predictor
#'
#' `step_lead_lag` creates a *specification* of a recipe step that
#'   will add new columns that are shifted forward (lag) or backward (lead).
#'   Data will by default include NA values where the shift was induced.
#'   These can be removed with [recipes::step_naomit()]. Samples should be ordered and
#'   have regular spacing (i.e. regular time series, regular spatial sampling).
#'
#' @inheritParams recipes::step_lag
#' @inheritParams recipes::step_center
#' @inheritParams step_distributed_lag
#' @param lag A vector of integers. Each specified column will be
#'  lagged for each value in the vector. Negative values are accepted and indicate
#'  leading the vector (i.e. the reverse of lagging)
#' @param n_subset A single integer. Subset every `n_subset` values.
#' @param n_shift A single integer amount to shift results in number of observations.
#' @param prefix A prefix for generated column names, default to "lag_lead_".
#' @param columns A character string of variable names that will
#'  be populated (eventually) by the `terms` argument.
#'
#' @details This step assumes that the data are already _in the proper sequential
#'  order_ for lagging. This `step` allows a vector to be shifted
#'  forward (lag) or backward (lead). While forward shifts are commonly used for
#'  lagged responses, there are cases where a backward shift may be useful. This
#'  can arise when there are unknown clock errors between two sensors making the
#'  response appear to occur before the input.  Another situation where a
#'  backward shift may be useful is in cyclical signals where alignment is
#'  unknown. The data can also efficiently be subsetted during the lag/leading
#'  process resulting in smaller model inputs while still utilizing the entire
#'  lag/lead history.
#'
#' @family row operation steps
#' @export
#' @rdname step_lead_lag
#'
#' @return An updated version of recipe with the new step added to the sequence
#'  of any existing operations.
#'
#' @examples
#' data(wipp30)
#'
#' recipe(wl~., data = wipp30) |>
#'   step_lead_lag(baro, lag = -2:2, n_subset = 1, n_shift = 0) |>
#'   prep()
#'
#' recipe(wl~ ., data = wipp30) |>
#'   step_lead_lag(baro, lag = -2:2, n_subset = 2, n_shift = 0) |>
#'   prep()
#'
#' recipe(wl~ ., data = wipp30) |>
#'   step_lead_lag(baro, lag = -2:2, n_subset = 2, n_shift = 1) |>
#'   prep()
#'
#' @seealso [recipes::step_lag()]
step_lead_lag <-
  function(recipe,
           ...,
           role = "predictor",
           trained = FALSE,
           lag = 1,
           n_subset = 1,
           n_shift = 0,
           prefix = "lead_lag_",
           keep_original_cols = FALSE,
           columns = NULL,
           skip = FALSE,
           id = rand_id("lead_lag")) {

    if(length(unique(lag)) < length(lag)) {
      rlang::warn("step_lead_lag should have uniquely valued 'lag'.  Taking unique values")
      lag <- unique(lag)
    }

    if(n_subset <= 0) {
      rlang::abort("'n_subset' in step_lead_lag should be greater than 0 in step_lead_lag")
    }

    if(n_shift >= n_subset) {
      rlang::abort("'n_shift' should be less than 'n_subset' in step_lead_lag")
    }


    add_step(
      recipe,
      step_lead_lag_new(
        terms = enquos(...),
        role = role,
        trained = trained,
        lag = lag,
        n_subset = n_subset,
        n_shift = n_shift,
        prefix = prefix,
        keep_original_cols = keep_original_cols,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_lead_lag_new <-
  function(terms, role, trained, lag, n_subset, n_shift, prefix, keep_original_cols, columns, skip, id) {
    step(
      subclass = "lead_lag",
      terms = terms,
      role = role,
      trained = trained,
      lag = lag,
      n_subset = n_subset,
      n_shift = n_shift,
      prefix = prefix,
      keep_original_cols = keep_original_cols,
      columns = columns,
      skip = skip,
      id = id
    )
  }

#' @export
prep.step_lead_lag <- function(x, training, info = NULL, ...) {

  col_names <- recipes_eval_select(x$terms, training, info)
  check_type(training[, col_names])

  step_lead_lag_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    lag = x$lag,
    n_subset = x$n_subset,
    n_shift = x$n_shift,
    prefix = x$prefix,
    keep_original_cols = x$keep_original_cols,
    columns = col_names,
    skip = x$skip,
    id = x$id
  )

}

#' @export
bake.step_lead_lag <- function(object, new_data, ...) {

  if(any(abs(object$lag) > nrow(new_data))) {
    rlang::abort("lag (lead) values cannot be greater than the number of rows in the dataset.")
  }

  if (!all(object$lag == as.integer(object$lag)))
    rlang::abort("step_lead_lag requires 'lead_lag' argument to be integer valued.")

  lag_mat <- lag_matrix(
    x = new_data[[object$columns]],
    lags = object$lag,
    n_subset = object$n_subset,
    n_shift = object$n_shift,
    var_name = object$prefix)
  colnames(lag_mat) <- paste0(colnames(lag_mat), '_', object$columns)


  if(object$n_subset > 1) {
    ind <- seq(object$n_shift + 1,
               nrow(new_data),
               object$n_subset)

    new_data <- bind_cols(new_data[ind,], lag_mat)
  } else {
    new_data <- bind_cols(new_data, lag_mat)
  }



  keep_original_cols <- get_keep_original_cols(object)
  if (!keep_original_cols) {
    new_data <-
      new_data[, !(colnames(new_data) %in% object$columns), drop = FALSE]
  }

  new_data
}

#' @importFrom recipes printer
print.step_lead_lag <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("lead_lag ",  sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }

#' @export
tidy.step_lead_lag <- function(x, ...) {
  tidy2.step_lead_lag(x, ...)
}

#' @rdname tidy2.recipe
#' @export
tidy2.step_lead_lag <- function(x, ...) {
  n_terms <- length(x$terms)
  res <-
    tibble(terms = rep(sel2char(x$terms), each = length(x$lag)),
           shift = rep(x$lag, times = n_terms)
    )
  res$key <- paste0(x$prefix, rep(x$lag, times = n_terms), '_', res$terms)
  res$id <- x$id
  res$step_name <- 'step_lead'
  res

}



