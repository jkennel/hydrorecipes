#' run_glmnet
#'
#' @param rec a prepped recipe (recipe)
#' @param new_data data to use with the recipe (data.frame)
#' @param outcome_column name of the outcome variable (character)
#' @param ...
#'
#' @return a list that contains the model fit info, the response weights,
#'  and the contributions for each component
#'
#' @export
#'
#' @examples
#' # kernel
#'
#' library(glmnet)
#' set.seed(1)
#' n_data   <- 200
#' kern_len <- 20
#' kern     <- 0.8^(0:kern_len) * 0.5
#' kern     <- kern / sum(kern)
#' input    <- rnorm(n_data)
#' outcome  <- fftw_convolve(input, rev(kern), align = 'right')
#'
#' plot(input, type = 'l', col = '#b47846', lwd = 2, ylab = 'values',
#'   main = 'orange is input, blue is output')
#' points(outcome, type = 'l', col = 'steelblue', lwd = 2)
#' data <- data.frame(outcome = head(outcome, 179), input = head(input, 179))
#' rec <- recipe(outcome~input, data) |>
#'   step_distributed_lag(input, knots = log_lags(4, 20)) |>
#'   step_rm('input') |>
#'   prep()
#' res <- run_glmnet(rec, new_data = data,
#'   outcome_column = 'outcome',
#'   alpha = 1,
#'   relax = TRUE)
#' plot(co~lag, res$response[[1]], type = 'o')
#' points(y = kern, x = 0:20, type = 'o', pch = 20, col = '#b47846', lwd = 2)
#' plot(outcome~predicted, res$decomp, type ='p', pch = 20, cex = 0.5)
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
  decomp <- as.data.frame(decomp)



  # calculate responses
  resp <- response(fit, rec)


  list(fit = fit,
       response = resp,
       decomp = decomp)


}
