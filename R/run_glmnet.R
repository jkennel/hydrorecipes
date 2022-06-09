#' run_glmnet
#'
#' @param rec
#' @param new_data
#' @param outcome_column
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
run_glmnet <- function(rec,
                       new_data = NULL,
                       outcome_column = 'outcome', ...) {

  # dataset for regression
  glm_rec <- rec |>
    recipes::bake(new_data = new_data)

  # save the complete cases
  complete <- complete.cases(glm_rec)

  # generate input matrices for glm_rec
  glm_in <- prep_for_glmnet(glm_rec,
                            outcome_column = outcome_column,
                            na_rm = FALSE)

  # fit model
  fit <- cv.glmnet(x = glm_in$x[complete, , drop = FALSE],
                   y = glm_in$y[complete, drop = FALSE],
                   ...)


  # calculate the contribution of each step
  decomp <- predict_terms(fit, rec, glm_in$x)
  decomp <- append(decomp, list(outcome = glm_in$y))
  decomp <- append(decomp, list(predicted = predict(fit, glm_in$x)))
  decomp <- data.table::as.data.table(decomp)



  # calculate responses
  resp <- response(fit, rec)


  list(fit = fit,
       response = resp,
       decomp = decomp)


}
