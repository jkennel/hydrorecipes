#' get_formula_vars
#'
#' @inheritParams lm
#'
#' @return
#' @export
#'
get_formula_vars <- function(formula, data) {

  left  <- rlang::f_lhs(formula)
  right <- rlang::f_rhs(formula)
  sym_dot <- sym(".")

  # check special cases
  if(left != sym_dot) {
    left  <- parse_formula(left)
    if(right == sym_dot) {
      right <- setdiff(names(data), left)
    } else {
      right <- parse_formula(right)
    }
  } else if (right == sym_dot) {
    right <- names(data)
    left  <- names(data)
  } else {
    right  <- parse_formula(right)
    left   <- setdiff(names(data), right)
  }


  list(
    predictors = setdiff(right, '+'),
    outcomes = setdiff(left, '+')
  )

}


parse_formula <- function(y){
  vapply(y,
         FUN = function(x) paste0(deparse(x)),
         FUN.VALUE = character(1L))
}



get_formula_vars_2 <- function(formula, data) {

  left  <- rlang::f_lhs(formula)
  right <- rlang::f_rhs(formula)
  sym_dot <- sym(".")

  nms <- names(data)

  # check special cases
  if(length(left) == 0L) {
    left <- ""
  }
  if(length(right) == 0L) {
    right <- ""
  }

  if(left != sym_dot) {
    left  <- parse_formula_2(left)

    if(right == sym_dot) {
      right <- setdiff(nms, left)
    } else {
      right <- parse_formula_2(right)
    }

  # both sides are "."
  } else if (right == sym_dot) {
    right <- nms
    left  <- nms
  # left side is "."
  } else {
    right  <- parse_formula_2(right)
    left   <- setdiff(nms, right)
  }


  list(
    predictors = intersect(nms, right),
    outcomes = intersect(nms, left)
  )

}

parse_formula_2 <- function(y){
  setdiff(unlist(strsplit(deparse(y), " +"), use.names = FALSE), "+")
}


# get the first class item
get_types <- function(data) {
  vapply(data,
         FUN = function(x) class(x)[1L],
         FUN.VALUE = character(1L))
}

# get column names
get_terms <- function(x) {
  vapply(x,
         FUN = rlang::as_name,
         FUN.VALUE = character(1L))
}


#' Make a random identification field for steps
#'
#' @export
#' @param prefix A single character string
#' @param len An integer for the number of random characters
#' @return A character string with the prefix and random letters separated by
#'  and underscore.
#'
#' @useDynLib frecipes, .registration = TRUE
#' @importFrom R6 R6Class
#' @importFrom Rcpp sourceCpp
#' @keywords internal
rand_id <- function(prefix = "step", len = 5L) {
  candidates <- c(letters, LETTERS, paste(0:9))
  paste(prefix,
        paste0(sample(candidates, len, replace = TRUE), collapse = ""),
        sep = "_"
  )
}
