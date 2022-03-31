#' prep_for_glmnet
#'
#' @param rec a recipe
#' @param outcome_column the column name that holds the outcome variable
#' @param na_rm remove incomplete rows
#'
#' @return
#' @export
#'
#' @examples
prep_for_glmnet <- function(rec,
                            outcome_column = 'outcome',
                            na_rm = TRUE) {

  r <- rec[, !names(rec) %in% 'predictor']
  y <- r[[outcome_column]]

  x <- r[, !names(r) %in% outcome_column]
  x <- do.call(cbind, x)

  if(na_rm) {
    complete <- complete.cases(r)
    y <- y[complete, drop = FALSE]
    x <- x[complete, , drop = FALSE]
  }

  list(x = x, y = y)

}
