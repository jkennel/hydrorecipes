#' @title tidy2
#'
#' @description  Turn an object into a tidy2 tibble
#'
#' @name tidy2

#' @param x An object to be converted into a tidy [tibble::tibble()].
#' @param ... Additional arguments to tidying method.
#' @return A [tibble::tibble()] with information about model components.
#'
#'
#' @export
tidy2 <- function(x, ...) {
  UseMethod("tidy2")
}

#' @title tidy2.recipe
#' @description `tidy2` will return a data frame that contains information
#'  regarding a recipe or operation within the recipe (when a `tidy`
#'  method for the operation exists).
#'
#' @name tidy2.recipe
#'
#'
#' @param x A `recipe` object, step, or check (trained or otherwise).
#' @param number An integer or `NA`. If missing and `id` is not provided,
#'  the return value is a list of the operations in the recipe.
#'  If a number is given, a `tidy2` method is executed for that operation
#'  in the recipe (if it exists). `number` must not be provided if
#'  `id` is.
#' @param id A character string or `NA`. If missing and `number` is not provided,
#'  the return value is a list of the operations in the recipe.
#'  If a character string is given, a `tidy2` method is executed for that
#'  operation in the recipe (if it exists). `id` must not be provided
#'  if `number` is.
#' @param ... Not currently used.
#' @return A tibble with columns that vary depending on what
#'  `tidy2` method is executed. When `number` and `id` are `NA`, a
#'  tibble with columns `number` (the operation iteration),
#'  `operation` (either "step" or "check"),
#'  `type` (the method, e.g. "nzv", "center"), a logical
#'  column called `trained` for whether the operation has been
#'  estimated using `prep`, a logical for `skip`, and a character column `id`.
NULL

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
    step_names <- vapply(x$steps, function(x) paste(x$columns, collapse = ''), character(1))

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
                  'step_lag',
                  'step_lead')) {
      res <- tidy2(x$steps[[number]], ...)
    }

  }

}

#' @export
tidy2.step <- function(x, ...) {
  rlang::abort(
    paste0(
      "No `tidy` method for a step with classes: ",
      paste0(class(x), collapse = ", ")
    )
  )
}

#' @export
tidy2.check <- function(x, ...) {
  rlang::abort(
    paste0(
      "No `tidy` method for a check with classes: ",
      paste0(class(x), collapse = ", ")
    )
  )
}





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
tidy2.step_intercept <- function(x, ...) {
  if (is_trained(x)) {
    terms <- names(x$objects)
  } else {
    terms <- sel2char(x$terms)
  }
  new_cols <- ncol(x$objects[[1]])

  ret <- tibble(terms = terms,
                id = x$id,
                step_name = 'step_intercept')
  ret$key <- paste(terms, "intercept", sep = "_")
  ret
}



