#' prep_for_glmnet
#'
#' @param rec
#' @param col_names
#'
#' @return
#' @export
#'
#' @examples
prep_for_glmnet <- function(rec, col_names = NULL) {

  r <- rec[, !names(rec) %in% 'predictor']
  r <- na.omit(r) # complete cases
  y <- r$outcome

  x <- r[, !names(r) %in% 'outcome']
  x <- do.call(cbind, x)

  list(x = x, y = y)

}
