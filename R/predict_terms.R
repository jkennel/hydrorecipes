#' predict_terms
#'
#' @description Predict the contribution for each step.
#'
#' @param fit A model object that has a `coefficients` method (e.g. `lm`)
#' @param rec A prepped `recipe`
#' @param data A data.frame with feature columns
#' @param ... Currently not used
#'
#' @return A data.frame of predicted values for each step (component contribution).
#'
#' @export
#' @examples
#' data(transducer)
#' transducer$datetime_num <- as.numeric(transducer$datetime)
#'
#' rec_toll_rasmussen <- recipe(wl ~ baro + et + datetime_num, transducer) |>
#'    step_lead_lag(baro, lag = log_lags(100, 86400 * 2 / 120)) |>
#'    step_ns(datetime_num, deg_free = 10) |>
#'    prep()
#'
#' input_toll_rasmussen <- rec_toll_rasmussen |> bake(new_data = NULL)
#'
#' fit_toll_rasmussen <- lm(wl ~ ., input_toll_rasmussen)
#' pred <- predict_terms(fit_toll_rasmussen,
#'                       rec_toll_rasmussen,
#'                       input_toll_rasmussen)
predict_terms <- function(fit, rec, data, ...) UseMethod("predict_terms")


#' @rdname predict_terms
#' @export
predict_terms.lm <- function(fit, rec, data, ...) {

  co <- coefficients(fit)

  predict_terms.numeric(co, rec, data, ...)

}


#' @rdname predict_terms
#' @export
#'
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
    } else {
      term_val <- rep(0, nrow(data))
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
