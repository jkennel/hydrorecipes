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

  ti <- ti[ti$variable %in% nms, ]

  x <- list()
  # create regression matrices
  x$term_info <- ti[ti$roles == id_type, ]
  x$term_info <- x$term_info[x$term_info$variable %in% unlist(vars), ]

  if (nrow(x$term_info) == 0) {
    stop(paste("Provided formula does not have any valid", id_type))
  }

  x$term_info$inds <- seq_len(nrow(x$term_info))

  x$term_info$ids <- which(nms %in% x$term_info$variable)

  x$to_rem <- collapse::missing_cases(new_data)
  x$data <- collapse::qM(unclass(new_data)[x$term_info$ids])
  x

}




# y = outcomes
# x = predictors
determine_coefficients <- function(x, y) {

  # solve
  fit <- llt_solve(
    x$data[!x$to_rem, , drop = FALSE],
    y$data[!y$to_rem, , drop = FALSE]
  )

  colnames(fit) <- y$term_info$variable
  rownames(fit) <- x$term_info$variable
  fit

}



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
predict_groups <- function(x, fit, steps) {

  # subsets are the regressor groups
  subsets <- subset_groups(x$term_info)
  lst <- list()

  for (i in seq_along(subsets)) {

    step_index <- unique(x$term_info[subsets[[i]], "step_index"])
    nms_vars   <- paste(steps[[step_index]]$columns, collapse = "_")
    step_name  <- unique(x$term_info[subsets[[i]], "step_name"])


    if (nms_vars == "" | step_name == "step_add_vars") {
      nms <- paste(colnames(fit), step_name, sep = "_")
    } else {
      nms <- paste(colnames(fit), step_name, nms_vars, sep = "_")
    }

    lst[[i]] <- collapse::mctl(
      x$data[, subsets[[i]], drop = FALSE] %*%
        fit[subsets[[i]], , drop = FALSE]
    )

    names(lst[[i]]) <- nms

  }

  lst
}

return_type <- function(x, type = "df") {

  # return types
  switch(
    type,
    "df" = collapse::qDF(x),
    "dt" = collapse::qDT(x),
    "tbl" = collapse::qTBL(x),
    "m" = collapse::qM(x),
    x
  )

}
