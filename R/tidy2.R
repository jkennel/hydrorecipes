#' @title tidy2
#'
#' @description  Turn an object into a tidy2 tibble.
#'
#' @name tidy2
#' @inheritParams recipes::tidy
#' @inheritParams generics::tidy
#'
#' @return A [tibble::tibble()] with information about model components.
#'
#'
#' @export
tidy2 <- function(x, ...) {
  UseMethod("tidy2")
}

#' @title tidy2.recipe
#' @description `tidy2` will return a data frame that contains information
#'  regarding a recipe or operation within the recipe (when a `tidy2`
#'  method for the operation exists). This method ensures that relevant data
#'  for `predict_terms` and `response` can be easily accessed from a recipe
#'  formulation.
#'
#' @name tidy2.recipe
#' @inheritParams recipes::tidy.recipe
#'
#' @return A tibble with columns that vary depending on what
#'  `tidy2` method is executed. When `number` and `id` are `NA`, a
#'  tibble with columns `number` (the operation iteration),
#'  `operation` (either "step" or "check"),
#'  `type` (the method, e.g. 'lead_lag', 'distributed_lag'), a logical
#'  column called `trained` for whether the operation has been
#'  estimated using `prep`, a logical for `skip`, and a character column `id`.
#'
NULL

#' rdname tidy2.recipe
#' @export
tidy2.recipe <- function(x, number = NA, id = NA, ...) {

  # add id = NA as default. If both ID & number are non-NA, error.
  # If number is NA and ID is not, select the step with the corresponding
  # ID. Only a single ID is allowed, as this follows the convention for number
  num_oper <- length(x$steps)
  pattern <- "(^step_)|(^check_)"

  if (!is.na(id)) {
    if (!is.na(number))
      rlang::abort("You may specify `number` or `id`, but not both.")
    if (length(id) != 1L && !is.character(id))
      rlang::abort("If `id` is provided, it must be a length 1 character vector.")
    step_ids <- vapply(x$steps, function(x) x$id, character(1))
    if(!(id %in% step_ids)) {
      rlang::abort("Supplied `id` not found in the recipe.")
    }
    number <- which(id == step_ids)
  }

  if (is.na(number)) {
    skipped <- vapply(x$steps, function(x) x$skip, logical(1))
    ids <- vapply(x$steps, function(x) x$id, character(1))
    step_names <- vapply(x$steps, function(x) paste(unique(tidy(x)$terms), collapse = ''), character(1))

    oper_classes <- lapply(x$steps, class)
    oper_classes <- grep("_", unlist(oper_classes), value = TRUE)

    oper <- strsplit(oper_classes, split = "_")
    oper <- vapply(oper, function(x) x[1], character(1))

    oper_types <- gsub(pattern, "", oper_classes)
    is_trained <- vapply(x$steps,
                         function(x) x$trained,
                         logical(1))
    res <- tibble(number = seq_along(x$steps),
                  operation = oper,
                  type = oper_types,
                  trained = is_trained,
                  skip = skipped,
                  id = ids,
                  step_name = step_names)
  } else {
    if (number > num_oper || length(number) > 1)
      rlang::abort(
        paste0(
          "`number` should be a single value between 1 and ",
          num_oper,
          "."
        )
      )

    res <- tidy(x$steps[[number]], ...)
    nm  <- class(x$steps[[number]])[1]

    if (nm %in% c('step_distributed_lag',
                  'step_earthtide',
                  'step_ns',
                  'step_intercept',
                  'step_lead_lag')) {

      res <- tidy2(x$steps[[number]], ...)

    }

  }
  res
}


#' @rdname tidy2.recipe
#' @export
tidy2.step <- function(x, ...) {
  rlang::abort(
    paste0(
      "No `tidy2` method for a step with classes: ",
      paste0(class(x), collapse = ", ")
    )
  )
}

#' @rdname tidy2.recipe
#' @export
tidy2.check <- function(x, ...) {
  rlang::abort(
    paste0(
      "No `tidy2` method for a check with classes: ",
      paste0(class(x), collapse = ", ")
    )
  )
}


#' @rdname tidy2.recipe
#' @export
tidy2.step_ns <- function(x, ...) {
  if (is_trained(x)) {
    terms <- names(x$objects)
  } else {
    terms <- sel2char(x$terms)
  }
  new_cols <- ncol(x$objects[[1]])

  ret <- tibble(terms = rep(terms, each = new_cols),
                id = x$id,
                step_name = 'step_ns')

  ret$key <- paste(rep(terms, each = new_cols), "ns",
                   rep(names0(new_cols, ""),
                       times = length(terms)), sep = "_")
  ret
}

#' @export
tidy.step_intercept <- function(x, ...) {
  tidy2.step_intercept(x, ...)
}

#' @rdname tidy2.recipe
#' @export
tidy2.step_intercept <- function(x, ...) {

  ret <- tibble(id = x$id,
                step_name = 'step_intercept')
  ret$key <- "intercept"
  ret
}



