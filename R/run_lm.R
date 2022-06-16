#' run_lm
#'
#' @param rec a prepped recipe (recipe)
#' @param new_data data to use with the recipe (data.frame)
#' @param ... arguments to pass to lm
#'
#' @return a list that contains the model fit info, the response weights,
#'  and the contributions for each component
#'
#' @export
#'
#' @examples
#' # kernel
#' library(hydrorecipes)
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
#' data <- data.frame(outcome = outcome, input = input)
#' rec <- recipe(outcome~input, data) |>
#'   step_distributed_lag(input, knots = log_lags(4, 20)) |>
#'   step_rm('input') |>
#'   prep()
#' res <- run_lm(rec,
#'   new_data = data)
#' plot(co~lag, res$response[[1]], type = 'o')
#' points(y = kern, x = 0:20, type = 'o', pch = 20, col = '#b47846', lwd = 2)
#' plot(outcome~predicted, res$components, type ='p', pch = 20, cex = 0.5)
run_lm <- function(rec,
                   new_data = NULL,
                   ...) {

  outcome_column <- rec$var_info
  outcome_column <- outcome_column$variable[
    which(outcome_column$role == 'outcome')]

  # dataset for regression
  lm_rec <- rec |>
    recipes::bake(new_data = new_data)


  # save the complete cases
  complete <- complete.cases(lm_rec)


  # fit model
  form <- formula(paste0(outcome_column, '~.'))
  fit  <- lm(formula = form, data = lm_rec, ...)


  # calculate the contribution of each step
  decomp <- predict_terms(fit, rec, lm_rec[complete,])
  decomp <- append(decomp, list(new_data[complete,]))
  decomp <- append(decomp, list(predicted = as.numeric(predict(fit, newx = lm_rec[complete,]))))
  decomp <- as.data.frame(decomp)



  # calculate responses
  resp <- response(fit, rec)


  resp <- list(fit = fit,
               response = resp,
               components = decomp)

  class(resp) <- 'deconvolution'
  resp

}
