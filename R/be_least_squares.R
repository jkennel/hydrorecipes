#' be_least_squares
#'
#' Calculate the barometric efficiency by using least squares.
#'
#' @inheritParams be_least_squares_diff
#'
#' @return barometric efficiency calculated by least squares
#' @export
#'
#'
#' @examples
#' library(data.table)
#' datetime <- seq.POSIXt(as.POSIXct("2016-01-01 12:00:00"),
#'                        as.POSIXct("2016-01-05 12:00:00"), by='hour' )
#' baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
#' wl <- -0.4 * baro
#' dat <- data.table(baro, wl, datetime)
#'
#' be_least_squares(dat)
#'
be_least_squares <- function(dat,
                             dep = 'wl',
                             ind = 'baro',
                             inverse = TRUE,
                             return_model = FALSE) {



  # don't modify existing data.table
  dat <- dat[, c(dep, ind), with = FALSE]


  if (inverse) {
    dat[, (dep) := -get(dep)]
  }

  fit_least_squares(dat, dep, ind, return_model)

}


#' fit_least_squares
#'
#' Least squares fit.
#'
#' @inheritParams be_least_squares_diff
#'
#' @return barometric efficiency calculated by least squares
#'
#'
fit_least_squares <- function(dat,
                   dep = 'wl',
                   ind = 'baro',
                   return_model = FALSE) {

  # fit regression
  frm <- formula(paste0(dep, "~", ind))
  be  <- lm(frm, dat)


  if (return_model) {
    return(be)
  }

  # calculate the slope
  return(as.numeric(coefficients(be)[2]))

}

#' be_least_squares_diff
#'
#' Calculate the barometric efficiency by using the least squares with differences
#'
#' @param dat data that has the independent and dependent variables (data.table)
#' @param dep name of the dependent variable column (character).  This is
#' typically the name for the column holding your water level data.
#' @param ind name of the independent variable column (character).  This is
#' typically the name for the column holding your barometric pressure data.
#' @param inverse  whether the barometric relationship is inverse
#' (TRUE means that when the barometric pressure goes up the measured water
#' level goes down (vented transducer, depth to water), FALSE means that when
#' the barometric pressure goes up so does the measured pressure
#' (non-vented transducer)) (logical).
#' @param lag_space space between difference calculation in number of observations
#' @param return_model whether to return the lm model or just the barometric/loading
#' efficiency (logical).#'
#' @return barometric efficiency calculated by least squares
#' @export
#'
#' @examples
#' library(data.table)
#' datetime <- seq.POSIXt(as.POSIXct("2016-01-01 12:00:00"),
#'                        as.POSIXct("2016-01-05 12:00:00"), by='hour' )
#' baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
#' wl <- -0.4 * baro
#' dat <- data.table(baro, wl, datetime)
#' be_least_squares_diff(dat, lag_space = 1)
#'
be_least_squares_diff <- function(dat,
                                  dep = 'wl',
                                  ind = 'baro',
                                  lag_space = 1,
                                  inverse = TRUE,
                                  return_model = FALSE) {


  # don't modify existing data.table
  dat <- dat[, c(dep, ind), with = FALSE]

  # calculate differences
  dat <- lag_difference(dat, dep, lag_space, inverse)
  dat <- lag_difference(dat, ind, lag_space)

  fit_least_squares(dat, dep, ind, return_model)

}

