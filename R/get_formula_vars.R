
# get_formula_vars_3 <- function(formula, data) {
#
#   left  <- rlang::f_lhs(formula)
#   right <- rlang::f_rhs(formula)
#   sym_dot <- as.symbol(".")
#
#   # check special cases
#   if (left != sym_dot) {
#     left  <- parse_formula(left)
#     if (right == sym_dot) {
#       right <- setdiff(names(data), left)
#     } else {
#       right <- parse_formula(right)
#     }
#   } else if (right == sym_dot) {
#     right <- names(data)
#     left  <- names(data)
#   } else {
#     right  <- parse_formula(right)
#     left   <- setdiff(names(data), right)
#   }
#
#
#   list(
#     predictors = setdiff(right, '+'),
#     outcomes = setdiff(left, '+')
#   )
#
# }
#
#
# parse_formula <- function(y){
#   vapply(y,
#          FUN = function(x) paste0(deparse(x)),
#          FUN.VALUE = character(1L))
# }
#
#
#
#
#
# get_formula_vars_2 <- function(formula, data) {
#
#   left  <- rlang::f_lhs(formula)
#   right <- rlang::f_rhs(formula)
#   sym_dot <- as.symbol(".")
#
#   nms <- names(data)
#
#   # check special cases
#   if (length(left) == 0L) {
#     left <- ""
#   }
#   if (length(right) == 0L) {
#     right <- ""
#   }
#
#   if (left != sym_dot) {
#     left  <- parse_formula_2(left)
#
#     if (right == sym_dot) {
#       right <- setdiff(nms, left)
#     } else {
#       right <- parse_formula_2(right)
#     }
#
#     # both sides are "."
#   } else if (right == sym_dot) {
#     right <- nms
#     left  <- nms
#     # left side is "."
#   } else {
#     right  <- parse_formula_2(right)
#     left   <- setdiff(nms, right)
#   }
#
#
#   list(
#     predictors = intersect(nms, right),
#     outcomes = intersect(nms, left)
#   )
#
# }


#' get_formula_vars
#'
#' @inheritParams lm
#'
#' @return
#' @export
#'
get_formula_vars <- function(formula, data) {

  dot <- "."

  form_char <- as.character(formula)

  left  <- all.vars(
    as.formula(file.path(form_char[2],
                         form_char[1], ".", fsep = ' ')), unique = FALSE)

  right <- all.vars(
    as.formula(file.path(form_char[3],
                         form_char[1], ".", fsep = ' ')), unique = FALSE)

  # remove the added "."
  left <- left[-length(left)]
  right <- right[-length(right)]

  nms <- names(data)

  if (any(right == dot) & any(left == dot)) {
    right <- nms
    left <- nms
  }
  if (any(right == dot) & !any(left == dot)) {
    right <- setdiff(nms, left)
  }
  if (!any(right == dot) & any(left == dot)) {
    left <- setdiff(nms, right)
  }

  list(
    predictors = intersect(nms, right),
    outcomes = intersect(nms, left)
  )

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# formula <- as.formula(x~.)
# data <- data.frame(x = 1, y = 3, z = 4, a = 1, b = 3)
# bench::mark(
#   get_formula_vars_3(formula, data),
#   get_formula_vars_2(formula, data)
# )


parse_formula_2 <- function(y){
  setdiff(unlist(strsplit(deparse(y), " +"), use.names = FALSE), "+")
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# get the first class item
get_sub_types <- function(data) {
  vapply(data,
         FUN = function(x) typeof(x),
         FUN.VALUE = character(1L),
         USE.NAMES = FALSE)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# get the first class item
get_types <- function(data) {
  vapply(data,
         FUN = function(x) {
           if (is.numeric(x)) {
             return("numeric")
           }
           if (is.character(x)) {
             return("character")
           }
           if (is.factor(x)) {
             return("factor")
           }
           if (is.logical(x)) {
             return("logical")
           }
           return("other")

         },
         FUN.VALUE = character(1L),
         USE.NAMES = FALSE)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# get column names
get_terms <- function(x) {
  vapply(x,
         FUN = rlang::as_name,
         FUN.VALUE = character(1L),
         USE.NAMES = FALSE)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



get_function_arguments <- function() {
  as.list(sys.frame(which = -1))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

get_function_arguments_no_rec <- function() {
  as.list(sys.frame(which = -1))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

get_terms_and_symbols <- function(terms) {
  # print(is.call(terms))
  # get function parameters to pass to parent
  lapply(terms, function(x) if (x != as.symbol("c")) x else NULL)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


get_terms_from_info <- function(terms, nms, info) {

  # for each call term, pull out the function selectors and pass the
  # necessary information
  include <- list()
  exclude <- list()

  # This loop needs to be refactored to something cleaner
  for (i in seq_along(terms)) {

    if (is.call(terms[[i]])) {
      # handle remove variable

      if (length(terms[[i]]) > 1 & as.character(terms[[i]])[1] == "-") {
        exclude[[i]] <- as.character(terms[[i]])[2]
      }  else {
        terms_list <- as.list(terms[[i]])

        # for contains
        if (length(terms_list) > 1) {
          include[[i]] <- (do.call(as.character(terms_list[[1]]),
                                   list(terms_list[[2]], nms)))
        } else {

          include[[i]] <- (do.call(as.character(terms[[i]]),
                                   list(nms, info)))
        }
      }

    } else {
      include[[i]] <- as.character(terms[[i]])
    }
  }
  # find matches for the data columns
  intersect(nms,
            setdiff(unique(unlist(include)),
                    unique(unlist(exclude))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


# type selectors
all_numeric <- function(nms, info) {
  nms[collapse::whichv(info$type, "numeric")]
}
non_numeric <- function(nms, info) {
  nms[collapse::whichv(info$type, "numeric", invert = FALSE)]
}
all_character <- function(nms, info) {
  nms[collapse::whichv(info$type, "character")]
}
all_datetime <- function(nms, info) {
  nms[collapse::whichv(info$type, "POSIXct")]
}

# sub type selectors
all_integer <- function(nms, info) {
  nms[collapse::whichv(info$sub_type, "integer")]
}
all_double <- function(nms, info) {
  nms[which(info$sub_type == "double")]
}
all_factor <- function(nms, info) {
  nms[which(info$sub_type == "factor")]
}
all_logical <- function(nms, info) {
  nms[which(info$sub_type == "logical")]
}
all_binary <- function(nms, info) {
  nms[collapse::whichv(info$sub_type, "binary")]
}
all_complex <- function(nms, info) {
  nms[collapse::whichv(info$sub_type, "complex")]
}

# role selectors
all_predictor <- function(nms, info) {
  nms[which(info$roles == "predictor")]
}
all_outcome <- function(nms, info) {
  nms[which(info$roles == "outcome")]
}



# grep selectors
contains <- function(to_find, nms) {
  nms[grepl(to_find, nms)]
}
not_contains <- function(to_find, nms) {
  nms[!grepl(to_find, nms)]
}
