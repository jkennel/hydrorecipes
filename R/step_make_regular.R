#' Create a lead predictor
#'
#' `step_make_regular` creates a *specification* of a recipe step that
#'   will generate an regularly spaced column and fill in the other columns with
#'   NA if missing.  It can also be used for subsetting.
#'
#' @inheritParams recipes::step_lag
#' @param delta A single numeric. This is the spacing for the regular series.
#' @param start A single numeric or value that can be converted to numeric. The
#' starting value of the regular series.  If left blank it defaults to the
#' minimum value.
#' @param prefix A prefix for generated column names, default to "lag_lead_".
#' @param columns A character string of variable names that will
#'  be populated (eventually) by the `terms` argument.
#'
#' @details This step sorts the data and converts a column to a regular series.
#' No interpolation is done and `NA` values are inserted where there is no
#' matching time.
#'
#' @family row operation steps
#' @export
#' @rdname step_make_regular
#'
#' @return An updated version of recipe with the new step added to the sequence
#'  of any existing operations.
#'
#' @examples
#' data(wipp30)
#' nrow(wipp30)
#' wipp30_missing <- wipp30[-(2:10),]
#' nrow(wipp30_missing)
#' wipp30_complete <- recipe(wl~., data = wipp30_missing) |>
#'   step_make_regular(time) |>
#'   prep() |>
#'   bake(new_data = NULL)
#' nrow(wipp30_complete)
#'
#' wipp30_subsample_2 <- recipe(wl~., data = wipp30_missing) |>
#'   step_make_regular(time, delta = 2) |>
#'   prep() |>
#'   bake(new_data = NULL)
#' nrow(wipp30_complete)
#'
#' @seealso [recipes::step_lag()] [step_distributed_lag()]
step_make_regular <-
  function(recipe,
           ...,
           role = "predictor",
           trained = FALSE,
           delta = 1,
           start = NULL,
           prefix = "make_regular_",
           columns = NULL,
           skip = FALSE,
           id = rand_id("make_regular")) {

    if(length(delta) != 1) {
      rlang::abort("step_make_regular requires a single 'delta' value")
    }


    add_step(
      recipe,
      step_make_regular_new(
        terms = enquos(...),
        role = role,
        trained = trained,
        delta = delta,
        start = start,
        prefix = prefix,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_make_regular_new <-
  function(terms, role, trained, delta, start, prefix, keep_original_cols, columns, skip, id) {
    step(
      subclass = "make_regular",
      terms = terms,
      role = role,
      trained = trained,
      delta = delta,
      start = start,
      prefix = prefix,
      columns = columns,
      skip = skip,
      id = id
    )
  }

#' @export
prep.step_make_regular <- function(x, training, info = NULL, ...) {

  col_names <- recipes_eval_select(x$terms, training, info)
  if(length(col_names) > 1) {
    rlang::abort("step_make_regular a single term")
  }
  check_type(training[, col_names])

  step_make_regular_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    delta = x$delta,
    start = x$start,
    prefix = x$prefix,
    columns = col_names,
    skip = x$skip,
    id = x$id
  )

}

#' @export
bake.step_make_regular <- function(object, new_data, ...) {

  # make sure data is sorted
  new_data <- new_data[order(new_data[[object$columns]]),]

  # create sequence
  rng <- as.numeric(range(new_data[[object$columns]], na.rm = TRUE))

  if(!is.null(object$start)) {
    if(object$start <= rng[2]){
      rng[1] <- object$start
    } else {
      rlang::abort('step_make_regular requires start to be less than the maximum value of the input column')
    }
  }
  new <- tibble(x = seq(rng[1], rng[2], by = object$delta))
  names(new) <- object$columns

  # do join
  left_join(new, new_data, by = object$columns)


}

#' @importFrom recipes printer
print.step_make_regular <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("make_regular ",  sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }

#' @export
tidy.step_make_regular <- function(x, ...) {
  tidy2.step_make_regular(x, ...)
}

#' @rdname tidy2.recipe
#' @export
tidy2.step_make_regular <- function(x, ...) {
  n_terms <- length(x$terms)
  res <-
    tibble(terms = sel2char(x$terms),
           delta = x$delta,
           start = x$start
    )
  res$key <- paste0(x$prefix, x$delta, '_', res$terms)
  res$id <- x$id
  res$step_name <- 'step_make_regular'
  res

}

#' @export
required_pkgs.step_make_regular <- function(x, ...) {
  c("hydrorecipes")
}
