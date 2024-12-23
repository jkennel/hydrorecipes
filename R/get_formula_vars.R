
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
#' @inheritParams stats::lm
#'
#' @return character vector of variable names
#' @export
#'
get_formula_vars <- function(formula, data) {

  dot <- "."
  if (is.null(formula)) {
    formula <- c("~", names(data)[1L], dot)
  }

  form_char <- as.character(formula)

  # print(form_char)
  left  <- all.vars(
    as.formula(file.path(form_char[2L],
                         form_char[1L], ".", fsep = ' ')), unique = FALSE)

  right <- all.vars(
    as.formula(file.path(form_char[3L],
                         form_char[1L], ".", fsep = ' ')), unique = FALSE)

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

select_fft_vars_list <- function(new_data, formula, columns) {

  vars_list <- names(new_data)

  if (!is.null(formula)) {
    vars_list <- get_formula_vars(formula = formula,
                                  data = unclass(new_data))

  } else {
    vars_list <- list(predictors = columns[-1L],
                      outcomes = columns[1L])
  }

  return(vars_list)
}


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
  as.list(sys.frame(which = -1L))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

get_function_arguments_no_rec <- function() {
  as.list(sys.frame(which = -1L))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

get_terms_and_symbols <- function(terms) {

  if(length(terms) > 1L) {
    if(terms[[1L]] == as.symbol("c")) {
      terms[[1L]] <- NULL
    }
  }

  terms
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


get_terms_from_info <- function(terms, nms) {

  # for each call term, pull out the function selectors and pass the
  # necessary information
  include <- list()
  exclude <- list()

  # This loop needs to be refactored to something cleaner
  for (i in seq_along(terms)) {

    if (is.call(terms[[i]])) {
      # handle remove variable

      if (length(terms[[i]]) > 1 & as.character(terms[[i]])[1L] == "-") {
        exclude[[i]] <- as.character(terms[[i]])[2L]
      }  else {
        terms_list <- as.list(terms[[i]])

        # for contains
        if (length(terms_list) > 1L) {
          include[[i]] <- (do.call(as.character(terms_list[[1L]]),
                                   list(terms_list[[2L]], nms)))
        } else {

          # include[[i]] <- (do.call(as.character(terms[[i]]),
          #                          list(nms, info)))
        }
      }

    } else {
      include[[i]] <- as.character(terms[[i]])
    }
  }

  # find matches for the data columns
  intersect(setdiff(collapse::funique(unlist(include)),
                    collapse::funique(unlist(exclude))),
            nms)

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
  nms[collapse::whichv(info$sub_type, "double")]
}
all_factor <- function(nms, info) {
  nms[collapse::whichv(info$sub_type, "factor")]
}
all_logical <- function(nms, info) {
  nms[collapse::whichv(info$sub_type, "logical")]
}
all_binary <- function(nms, info) {
  nms[collapse::whichv(info$sub_type, "binary")]
}
all_complex <- function(nms, info) {
  nms[collapse::whichv(info$sub_type, "complex")]
}

# role selectors
all_predictor <- function(nms, info) {
  nms[collapse::whichv(info$roles, "predictor")]
}
all_outcome <- function(nms, info) {
  nms[collapse::whichv(info$roles, "outcome")]
}



# grep selectors
contains <- function(to_find, nms) {
  nms[grepl(to_find, nms)]
}
not_contains <- function(to_find, nms) {
  nms[!grepl(to_find, nms)]
}
