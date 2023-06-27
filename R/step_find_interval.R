#' Create a lead predictor
#'
#' `step_find_interval` creates a *specification* of a recipe step that
#'   will add binary terms (one-hot) that are based on provided boundary values.
#'   This is an efficient combination of `recipes::step_cut` and
#'   `recipes::step_dummy`.
#'
#' @inheritParams recipes::step_dummy
#' @inheritParams recipes::step_cut
#' @inheritParams recipes::step_lag
#' @inheritParams base::findInterval
#' @param encoding Character string specifying the type of encoding to use. The
#'   options are "dummy", "one_hot", "factor" or "integer" encoding. Defaults to
#'   "one_hot".
#' @param prefix A prefix for generated column names, default to
#'   "find_interval_".
#'
#' @details This function uses `base::find_interval()` to find data within a set
#'  of provided intervals. The data can either be returned as an factor column
#'  or with dummy encoding.
#'
#' @export
#' @rdname step_find_interval
#'
#' @return An updated version of recipe with the new step added to the sequence
#'  of any existing operations.
#'
#' @examples
#' data(wipp30)
#'
#' r_fi_one_hot <- recipe(wl ~ ., data = wipp30) |>
#'   step_find_interval(time,
#'     vec = c(9000, 12000),
#'     keep_original_cols = TRUE,
#'     encoding = "one_hot"
#'   ) |>
#'   prep() |>
#'   bake(new_data = NULL)
#'
#' r_fi_dummy <- recipe(wl ~ ., data = wipp30) |>
#'   step_find_interval(time,
#'     vec = c(9000, 12000),
#'     encoding = "dummy"
#'   ) |>
#'   prep() |>
#'   bake(new_data = NULL)
#'
#' r_fact <- recipe(wl ~ ., data = wipp30) |>
#'   step_find_interval(time,
#'     vec = c(9000, 12000),
#'     encoding = "factor"
#'   ) |>
#'   prep() |>
#'   bake(new_data = NULL)
#'
#' @seealso [recipes::step_cut()] [recipes::step_dummy]
step_find_interval <-
  function(recipe,
           ...,
           role = "predictor",
           trained = FALSE,
           vec = NULL,
           encoding = "one_hot",
           prefix = "find_interval_",
           keep_original_cols = FALSE,
           columns = NULL,
           skip = FALSE,
           id = rand_id("find_interval")) {
    if (is.null(vec)) {
      rlang::abort("'vec' in step_find_interval cannot be NULL")
    }


    add_step(
      recipe,
      step_find_interval_new(
        terms = enquos(...),
        role = role,
        trained = trained,
        vec = vec,
        encoding = encoding,
        prefix = prefix,
        keep_original_cols = keep_original_cols,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_find_interval_new <-
  function(
      terms, role, trained, vec, encoding, prefix, keep_original_cols,
      columns, skip, id) {
    step(
      subclass = "find_interval",
      terms = terms,
      role = role,
      trained = trained,
      vec = vec,
      encoding = encoding,
      prefix = prefix,
      keep_original_cols = keep_original_cols,
      columns = columns,
      skip = skip,
      id = id
    )
  }

#' @export
prep.step_find_interval <- function(x, training, info = NULL, ...) {
  col_names <- recipes_eval_select(x$terms, training, info)
  check_type(training[, col_names])
  vec <- sort(x$vec)

  step_find_interval_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    vec = vec,
    encoding = x$encoding,
    prefix = x$prefix,
    keep_original_cols = x$keep_original_cols,
    columns = col_names,
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_find_interval <- function(object, new_data, ...) {
  for (i in seq_along(object$columns)) {
    column_name <- object$columns[i]
    interval <- findInterval(new_data[[column_name]], vec = object$vec)
    interval_f <- collapse::qF(interval)
    n_vec <- nlevels(interval_f)

    if (object$encoding %in% c("one_hot", "dummy")) {
      column_names <- paste0(object$prefix, column_name, "_", 0L:(n_vec - 1L))
    } else {
      column_names <- paste0(object$prefix, column_name)
    }

    if (n_vec <= 1) {
      if (object$encoding != "dummy") {
        new_data[, column_names] <- 1L
      }
    } else if (object$encoding == "one_hot") {
      interval_dummy <- to_dummy(interval_f, n_vec)
      colnames(interval_dummy) <- column_names
      new_data <- bind_cols(new_data, as_tibble(interval_dummy))
    } else if (object$encoding == "dummy") {
      interval_dummy <- to_dummy(interval_f, n_vec)
      colnames(interval_dummy) <- column_names
      new_data <- bind_cols(
        new_data, as_tibble(interval_dummy[, -1L, drop = FALSE])
      )
    } else if (object$encoding == "factor") {
      new_data[, column_names] <- interval_f
    } else {
      new_data[, column_names] <- interval
    }
  }



  keep_original_cols <- get_keep_original_cols(object)
  if (!keep_original_cols) {
    new_data <-
      new_data[, !(colnames(new_data) %in% object$columns), drop = FALSE]
  }

  new_data
}

#' @importFrom recipes printer
print.step_find_interval <-
  function(x, width = max(20L, options()$width - 30L), ...) {
    cat("find_interval ", sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }

#' @export
tidy.step_find_interval <- function(x, ...) {
  tidy2.step_find_interval(x, ...)
}

#' @rdname tidy2.recipe
#' @export
tidy2.step_find_interval <- function(x, ...) {
  n_terms <- length(x$terms)
  if (is_trained(x)) {
    term_names <- x$columns
  } else {
    term_names <- sel2char(x$terms)
  }
  res <-
    tibble(
      terms = rep(sel2char(x$terms), each = length(x$vec)),
      vec = rep(x$vec, times = n_terms),
      key = paste0(x$prefix, term_names)
    )
  res$id <- x$id
  res$step_name <- "step_find_interval"
  res
}

#' @export
required_pkgs.step_find_interval <- function(x, ...) {
  c("hydrorecipes", "collapse")
}
