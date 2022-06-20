#' predict_terms
#'
#' @description Predict the contribution for each step.
#'
#' @param fit a model object that has a `coefficients` method
#' @param rec a prepped `recipe`
#' @param data data.frame with feature columns
#' @param ... not used
#'
#' @return a list of response functions corresponding to each step
#' @export
#'
predict_terms <- function(fit, rec, data, ...) UseMethod("predict_terms")


#' @rdname predict_terms
#' @export
predict_terms.lm <- function(fit, rec, data, ...) {

  co <- coefficients(fit)

  predict_terms.numeric(co, rec, data, ...)

}


#' @rdname predict_terms
#' @export
predict_terms.cv.glmnet <- function(fit, rec, data, ...) {

  co <- coefficients(fit)
  co_names <- rownames(co)
  co <- as.vector(co)
  names(co) <- co_names

  predict_terms.numeric(co, rec, data, ...)

}


#' @rdname predict_terms
#' @export
predict_terms.numeric <- function(fit, rec, data, ...) {

  rec_steps <- tidy2(rec)
  rec_steps <- rec_steps[!rec_steps$type %in% 'rm', ]

  # output list
  resp <- vector(mode = "list",
                 length = nrow(rec_steps))
  names(resp) <- paste(rec_steps$type, rec_steps$step_name, sep = '_')


  for (i in 1:nrow(rec_steps)) {

    step_info <- tidy2(rec, i)

    co_sub <- get_coefficients(fit, step_info)

    if (length(co_sub) >= 1) {
      term_val <- as.matrix(data[, names(co_sub)]) %*% as.matrix(co_sub)
    }

    resp[[i]] <- term_val
  }

  df <- as.data.frame(resp)

  if('(Intercept)' %in% names(fit)) {
    df[['intercept']] <- fit['(Intercept)']
  }

  df$predicted <- rowSums(df)
  df

}
