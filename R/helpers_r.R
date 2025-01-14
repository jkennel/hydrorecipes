#' pad_num
#'
#' This function creates a sequence of numbers with 0 padding based on the
#' series length.
#'
#' @param n the number of columns
#' @param pad the character to use for the padding.
#'
#' @return a character string padded by "0"
#' @export
#'
pad_num <- function(n, pad = "0") {
  width <- floor(log10(n)) + 1L

  formatC(seq_len(n),
    width = width,
    format = "d",
    flag = "0"
  )
}

name_columns <- function(id, column_name, n) {

  if (is.null(column_name)) {
    if (n < 2) {
      return(file.path(id, fsep = "_"))
    }
    return(file.path(id, pad_num(n), fsep = "_"))
  }

  if (n < 2) {
    return(file.path(id, column_name, fsep = "_"))
  }


  file.path(id, column_name, pad_num(n), fsep = "_")
}




#' Make a random identification field for steps
#'
#' @export
#' @param prefix A single character string
#' @param len An integer for the number of random characters
#' @return A character string with the prefix and random letters separated by
#'  and underscore.
#'
#' @keywords internal
rand_id <- function(prefix = "step", len = 5L) {
  candidates <- c(letters, LETTERS, paste(0:9))
  paste(prefix,
    paste0(sample(candidates, len, replace = TRUE), collapse = ""),
    sep = "_"
  )
}

# parse_variables <- function(env, arg_nms) {
#
#       env_nms <- names(env)
#       super_nms <- arg_nms
#       env_nms_sub <- intersect(super_nms, env_nms)[-1L]
#
#       inputs <- c(
#         as.list(substitute(test, environment())),
#         as.list(environment())[env_nms_sub]
#       )
#
# }



# regression helpers ------------------------------------------------------
# predictors outcomes
get_regression_data <- function(new_data,
                                term_info,
                                vars,
                                id_type = "predictor") {

  nms <- unique(names(new_data))

  # term info data
  ti <- collapse::qDF(term_info)
  ti <- ti[ti$source != "removed", ]

  ti <- ti[ti$variable %iin% nms, ]

  x <- list()
  # create regression matrices
  x$term_info <- ti[ti$roles == id_type, ]
  x$term_info <- x$term_info[x$term_info$variable %iin% unlist(vars), ]

  if (nrow(x$term_info) == 0) {
    stop(paste("Provided formula does not have any valid", id_type))
  }

  x$term_info$inds <- seq_len(nrow(x$term_info))

  x$term_info$ids <- nms %iin% x$term_info$variable

  x$to_rem <- collapse::missing_cases(new_data)
  x$data <- collapse::qM(unclass(new_data)[x$term_info$ids])
  x

}



# y = outcomes
# x = predictors
determine_coefficients <- function(x, y, has_na, decomp) {

  # print(str(x[!has_na, , drop = FALSE]))
  # print(str(y[!has_na, , drop = FALSE]))
  # solve
  fit <- llt_solve_full(
    x[!has_na, , drop = FALSE],
    y[!has_na, , drop = FALSE],
    decomp
  )

  # f <- lm.fit(x[!has_na, , drop = FALSE],
  #             y[!has_na, , drop = FALSE])
  # print(fit$coefficients)
  # print(f$coefficients)
  # plot(f$fitted.values, type='l', lwd = 3)
  # points(fit$fitted.values, type = 'l', col = "green")
  colnames(fit$coefficients) <- colnames(y)
  rownames(fit$coefficients) <- colnames(x)
  fit

}

# X <- matrix(rnorm(1e7), ncol = 10)
# y <- as.matrix(rnorm(1e6))
# bench::mark(hydrorecipes:::llt_solve(X, y),
#             hydrorecipes:::llt_solve_full(X,y, list(a = c(0,2), b = c(2,2))),
#             lm.fit(X, y), check = FALSE)

subset_groups <- function(x) {
  split(
    x$inds,
    data.table::rleid(x$step_index)
  )
}

response_groups <- function(steps, x, fit) {
  # subsets are the regressor groups
  subsets <- subset_groups(x$term_info)

  lst <- list()
  for (i in seq_along(subsets)) {
    lst[[i]] <- steps[[i]]$response(fit[subsets[[i]], , drop = FALSE])
  }

  lst
}

# x = predictors
predict_groups <- function(x, fit, step_vars, step_names) {

  # subsets are the regressor groups
  lst <- list()

  for (i in seq_along(step_vars)) {

    nms_vars   <- paste(step_vars[[i]], collapse = "_")


    lst[[i]] <- collapse::mctl(
      x$data[, step_vars[[i]], drop = FALSE] %*%
        fit[step_vars[[i]], , drop = FALSE]
    )

    names(lst[[i]]) <- step_names

  }

  lst
}

predict_each_step <- function(x, fit, step_vars, step_names) {

  # subsets are the regressor groups
  lst <- list()

  predicted <- rep.int(0.0, nrow(x))

  for (i in seq_along(step_vars)) {

    nms <- step_vars[[i]]

    if (is.null(nms)) {
      next
    }

    wh <- colnames(x) %iin% nms

    if (!any(wh)) {
      next
    }


    nms_step  <- paste(paste(step_names[[i]], collapse = "_"),
                       paste(nms, collapse = "_"),
                       sep = "_")


    lst[i] <- collapse::mctl(x[, wh, drop = FALSE] %*% fit[wh, , drop = FALSE])
    predicted %+=% lst[[i]]

    names(lst)[i] <- nms_step
  }

  lst[["predicted"]] <- predicted
  lst
}


# formula can be used to subset or separate predictors and outcomes
return_type <- function(x, type = "df", formula = NULL, combined = TRUE) {

  if (!is.null(formula)) {
    vars_list <- get_formula_vars(formula = formula, data = unclass(x))
  } else {
    vars_list <- names(x)

    if (!combined) {
      vars_list <- list(predictors = vars_list[-1L],
                        outcomes = vars_list[1L])
    }

  }

  if (combined) {
    vars_list <- unlist(vars_list)

    # return types
    x <- switch(
      type,
      "df" = collapse::qDF(unclass(x)[vars_list]),
      "dt" = collapse::qDT(unclass(x)[vars_list]),
      "tbl" = collapse::qTBL(unclass(x)[vars_list]),
      "m" = collapse::qM(unclass(x)[vars_list]),
      unclass(x)[vars_list]
    )

    return(x)
  }

  x <- switch(
    type,
    "df" = list(predictors  = collapse::qDF(unclass(x)[vars_list[[1L]]]),
                outcomes    = collapse::qDF(unclass(x)[vars_list[[2L]]])),
    "dt" = list(predictors  = collapse::qDT(unclass(x)[vars_list[[1L]]]),
                outcomes    = collapse::qDT(unclass(x)[vars_list[[2L]]])),
    "tbl" = list(predictors = collapse::qTBL(unclass(x)[vars_list[[1L]]]),
                 outcomes   = collapse::qTBL(unclass(x)[vars_list[[2L]]])),
    "m"   = list(predictors = collapse::qM(unclass(x)[vars_list[[1L]]]),
                 outcomes   = collapse::qM(unclass(x)[vars_list[[2L]]])),

    list(predictors = unclass(x)[vars_list[[1L]]],
         outcomes   = unclass(x)[vars_list[[2L]]])

  )




}
