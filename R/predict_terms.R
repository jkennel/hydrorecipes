#' predict_terms
#'
#' @param fit a model object that has a `coefficients` method
#' @param rec a prepped recipe
#' @param data dataset with feature columns
#' @param ... other arguments
#'
#' @return a list of response functions corresponding to each step
#' @export
#'
predict_terms <- function(fit, rec, data, ...) UseMethod("predict_terms")


#' @rdname predict_terms
#' @export
predict_terms.lm <- function(fit, rec, data, ...) {

  co <- coefficients(fit)

  rec_steps <- tidy2(rec)
  rec_steps <- rec_steps[rec_steps$type != 'rm', ]

  resp <- vector(mode = "list",
                 length = nrow(rec_steps))
  names(resp) <- paste(rec_steps$type, rec_steps$step_name, sep = '_')

  for (i in 1:nrow(rec_steps)) {

    step_info <- tidy2(rec, i)

    type <- rec_steps$type[i]

    co_sub <- get_coefficients(co, step_info)

    if (length(co_sub) >= 1) {
      term_val <- as.matrix(data[, names(co_sub)]) %*% co_sub
      # plot(term_val, type = 'l')
    }
    resp[[i]] <- term_val
  }
  resp

}


#' @rdname predict_terms
#' @export
predict_terms.cv.glmnet <- function(fit, rec, data, ...) {

  co <- coefficients(fit)
  co_names <- rownames(co)
  co <- as.vector(co)
  names(co) <- co_names

  rec_steps <- tidy2(rec)
  rec_steps <- rec_steps[rec_steps$type != 'rm', ]

  # output list
  resp <- vector(mode = "list",
                 length = nrow(rec_steps))
  names(resp) <- paste(rec_steps$type, rec_steps$step_name, sep = '_')


  for (i in 1:nrow(rec_steps)) {

    step_info <- tidy2(rec, i)

    co_sub <- get_coefficients(co, step_info)

    if (length(co_sub) >= 1) {
      term_val <- as.matrix(data[, names(co_sub)]) %*% as.matrix(co_sub)
      # plot(term_val, type = 'l', main = rec_steps$id[i])
    }
    resp[[i]] <- term_val
  }
  resp

}

