#' Create a lead predictor
#'
#' `step_lead` creates a *specification* of a recipe step that
#'   will add new columns of lead data. lead data will
#'   by default include NA values where the lead was induced.
#'   These can be removed with [step_naomit()], or you may
#'   specify an alternative filler value with the `default`
#'   argument.
#'
#' @inheritParams recipes::step_lag
#' @inheritParams recipes::step_center
#' @param lead A vector of positive integers. Each specified column will be
#'  lead for each value in the vector.
#' @param prefix A prefix for generated column names, default to "lead_".
#' @param columns A character string of variable names that will
#'  be populated (eventually) by the `terms` argument.
#' @param default Passed to `dplyr::lead`, determines what fills empty rows
#'   left by leading (defaults to NA).
#' @details The step assumes that the data are already _in the proper sequential
#'  order_ for leading.
#' @family row operation steps
#' @export
#' @rdname step_lead
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
#'   step_lead(index, day, lead = 2:3) |>
#'   prep() |>
#'   bake(df)
step_lead <-
  function(recipe,
           ...,
           role = "predictor",
           trained = FALSE,
           lead = 1,
           prefix = "lead_",
           default = NA,
           columns = NULL,
           skip = FALSE,
           id = rand_id("lead")) {
    add_step(
      recipe,
      step_lead_new(
        terms = enquos(...),
        role = role,
        trained = trained,
        lead = lead,
        default = default,
        prefix = prefix,
        columns = columns,
        skip = skip,
        id = id
      )
    )
  }

step_lead_new <-
  function(terms, role, trained, lead, default, prefix, columns, skip, id) {
    step(
      subclass = "lead",
      terms = terms,
      role = role,
      trained = trained,
      lead = lead,
      default = default,
      prefix = prefix,
      columns = columns,
      skip = skip,
      id = id
    )
  }

#' @export
prep.step_lead <- function(x, training, info = NULL, ...) {
  step_lead_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    lead = x$lead,
    default = x$default,
    prefix = x$prefix,
    columns = recipes_eval_select(x$terms, training, info),
    skip = x$skip,
    id = x$id
  )
}

#' @export
bake.step_lead <- function(object, new_data, ...) {

  if (!all(object$lead == as.integer(object$lead)))
    rlang::abort("step_lead requires 'lead' argument to be integer valued.")

  make_call <- function(col, lead_val) {
    call2(
      "lead",
      x = sym(col),
      n = lead_val,
      default = object$default,
      .ns = "dplyr"
    )
  }

  grid <- tidyr::expand_grid(col = object$columns, lead_val = object$lead)
  calls <- purrr::map2(grid$col, grid$lead_val, make_call)
  newname <- as.character(glue::glue("{object$prefix}{grid$lead_val}_{grid$col}"))
  calls <- check_name(calls, new_data, object, newname, TRUE)

  as_tibble(mutate(new_data, !!!calls))
}

print.step_lead <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("Leading ",  sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }

#' @rdname tidy.recipe
#' @param x A `step_lead` object.
#' @export
tidy.step_lead <- function(x, ...) {
  n_terms <- length(x$terms)
  res <-
    tibble(terms = rep(sel2char(x$terms), each = length(x$lead)),
           shift = rep(-x$lead, times = n_terms)
    )
  res$key <- paste0(x$prefix, rep(x$lead, times = n_terms), '_', res$terms)
  res$id <- x$id
  res
}



