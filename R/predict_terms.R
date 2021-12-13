#' predict_terms
#'
#' @param fit a model object that has a `coefficients` method
#' @param rec a prepped recipe
#' @param data dataset with feature columns
#'
#' @return a list of response functions corresponding to each step
#' @export
#'
predict_terms <- function(fit, rec, data) {

  co <- coefficients(fit)

  rec_steps <- tidy(rec)

  resp <- vector(mode = "list", length = nrow(rec_steps))
  names(resp) <- rec_steps$type

  for (i in 1:nrow(rec_steps)) {

    step_info <- tidy(rec, i)

    type <- rec_steps$type[i]

    co_sub <- get_coefficients(co, step_info)

    if (length(co_sub) >= 1) {
      term_val <- as.matrix(data[, names(co_sub)]) %*% co_sub
      plot(term_val, type = 'l')
    }

  }

}
