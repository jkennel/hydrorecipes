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
  sym_dot <- as.symbol(".")

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
  sym_dot <- as.symbol(".")

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



get_formula_vars_3 <- function(formula, data) {

  dot <- "."

  form_char <- as.character(formula)

  left  <- all.vars(as.formula(file.path(form_char[2],
                                         form_char[1], ".", fsep = ' ')), unique = FALSE)
  right <- all.vars(as.formula(file.path(form_char[3],
                                         form_char[1], ".", fsep = ' ')), unique = FALSE)

  # remove the added "."
  left <- left[-length(left)]
  right <- right[-length(right)]

  nms <- names(data)

  if(any(right == dot) & any(left == dot)) {
    right <- nms
    left <- nms
  }
  if(any(right == dot) & !any(left == dot)) {
    right <- setdiff(nms, left)
  }
  if(!any(right == dot) & any(left == dot)) {
    left <- setdiff(nms, right)
  }

  list(
    predictors = intersect(nms, right),
    outcomes = intersect(nms, left)
  )

}

formula <- as.formula(x~.)
data <- data.frame(x = 1, y = 3, z = 4, a = 1, b = 3)
bench::mark(
  get_formula_vars_3(formula, data),
  get_formula_vars_2(formula, data)
)


parse_formula_2 <- function(y){
  setdiff(unlist(strsplit(deparse(y), " +"), use.names = FALSE), "+")
}


# get the first class item
get_types <- function(data) {
  vapply(data,
         FUN = function(x) class(x)[1L],
         FUN.VALUE = character(1L),
         USE.NAMES = FALSE)
}

# get column names
get_terms <- function(x) {
  vapply(x,
         FUN = rlang::as_name,
         FUN.VALUE = character(1L),
         USE.NAMES = FALSE)
}


