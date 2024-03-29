% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict_terms.R
\name{predict_terms}
\alias{predict_terms}
\alias{predict_terms.lm}
\alias{predict_terms.cv.glmnet}
\alias{predict_terms.numeric}
\title{predict_terms}
\usage{
predict_terms(fit, rec, data, ...)

\method{predict_terms}{lm}(fit, rec, data, ...)

\method{predict_terms}{cv.glmnet}(fit, rec, data, ...)

\method{predict_terms}{numeric}(fit, rec, data, ...)
}
\arguments{
\item{fit}{A model object that has a \code{coefficients} method (e.g. \code{lm})}

\item{rec}{A prepped \code{recipe}}

\item{data}{A data.frame with feature columns}

\item{...}{Currently not used}
}
\value{
A data.frame of predicted values for each step (component contribution).
}
\description{
Predict the contribution for each step.
}
\examples{
data(transducer)
transducer$datetime_num <- as.numeric(transducer$datetime)

rec_toll_rasmussen <- recipe(wl ~ baro + et + datetime_num, transducer) |>
   step_lead_lag(baro, lag = log_lags(100, 86400 * 2 / 120)) |>
   step_ns(datetime_num, deg_free = 10) |>
   prep()

input_toll_rasmussen <- rec_toll_rasmussen |> bake(new_data = NULL)

fit_toll_rasmussen <- lm(wl ~ ., input_toll_rasmussen)
pred <- predict_terms(fit_toll_rasmussen,
                      rec_toll_rasmussen,
                      input_toll_rasmussen)
}
