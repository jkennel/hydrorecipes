#' Create a lead predictor
#'
#' `step_lead_lag` creates a *specification* of a recipe step that
#'   will add new columns of lead data. lead data will
#'   by default include NA values where the lead was induced.
#'   These can be removed with [step_naomit()], or you may
#'   specify an alternative filler value with the `default`
#'   argument.
#'
#' @inheritParams recipes::step_lag
#' @inheritParams recipes::step_center
#' @param lag A vector of integers. Each specified column will be
#'  lag for each value in the vector. Negative values are accepted and indicate
#'  leading the vector (i.e. the reverse of lagging)
#' @param prefix A prefix for generated column names, default to "lag_lead_".
#' @param columns A character string of variable names that will
#'  be populated (eventually) by the `terms` argument.
#' @details The step assumes that the data are already _in the proper sequential
#'  order_ for leading.
#' @family row operation steps
#' @export
#' @rdname step_lead_lag
#'
#' @examples
#' n <- 10
#' start <- as.Date('1999/01/01')
#' end <- as.Date('1999/01/10')
#'
#' df <- data.frame(x = runif(n),
#'                  index = 1:n,
#'                  day = seq(start, end, by = "day"))
#'
#' recipe(~ ., data = df) |>
#'   step_lead_lag(index, lag = 1, n_subset = 1, n_shift = 0) |>
#'   prep() |>
#'   bake(df)
#'
step_lead_lag <-
  function(recipe,
           ...,
           role = "predictor",
           trained = FALSE,
           lag = 1,
           n_subset = 1,
           n_shift = 0,
           prefix = "lead_lag_",
           columns = NULL,
           skip = FALSE,
           id = rand_id("lead_lag")) {
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
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_lead_lag_new <-
  function(terms, role, trained, lag, n_subset, n_shift, prefix, columns, skip, id) {
    step(
      subclass = "lead_lag",
      terms = terms,
      role = role,
      trained = trained,
      lag = lag,
      n_subset = n_subset,
      n_shift = n_shift,
      prefix = prefix,
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
    columns = col_names,
    skip = x$skip,
    id = x$id
  )

}

#' @export
bake.step_lead_lag <- function(object, new_data, ...) {

  if (!all(object$lag == as.integer(object$lag)))
    rlang::abort("step_lead_lag requires 'lead_lag' argument to be integer valued.")


  if(object$n_subset > 1) {
    return(
      bind_cols(new_data[seq(object$n_shift+1, nrow(new_data) - object$n_shift, object$n_subset),],
                lag_matrix(
                  new_data[[object$columns]],
                  object$lag,
                  object$n_subset,
                  object$n_shift))
    )
  }


  bind_cols(new_data,
            lag_matrix(
              new_data[[object$columns]],
              object$lag,
              object$n_subset,
              object$n_shift))



}



#' @export
tidy.step_lead_lag <- function(x, ...) {
  tidy2(x, ...)
}


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



