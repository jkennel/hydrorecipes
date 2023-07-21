#' @title
#' be_clark
#'
#' @description
#' Clark 1967 solution for calculating barometric efficiency.
#'
#' @inheritParams be_least_squares_diff
#'
#' @return \code{lm} linear model for Clark's method.  The coefficient is the BE/LE
#'
#' @references Clark, W., 1967. Computing the Barometric Efficiency of a well. Proc. Am. Soc. Civ. Eng. J. Hydraul. Div. V. 93, 93–98.
#'
#' @import data.table
#' @export
#'
#' @examples
#' library(data.table)
#' datetime <- seq.POSIXt(as.POSIXct("2016-01-01 12:00:00"),
#'                        as.POSIXct("2016-01-05 12:00:00"), by='sec' )
#' baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
#' wl   <- -0.4 * baro + rnorm(length(datetime), sd = 0.02)
#' dat  <- data.table(baro, wl, datetime)
#' be_clark(dat, dep = 'wl', ind = 'baro', lag_space = 1, inverse = TRUE)
#'
be_clark <- function(dat,
                     dep = 'wl',
                     ind = 'baro',
                     lag_space = 1,
                     inverse = TRUE,
                     return_model = FALSE) {


  # calculate differences
  dat <- dat[, c(dep, ind), with = FALSE]
  dat <- lag_difference(dat, dep, lag_space, inverse)
  dat <- lag_difference(dat, ind, lag_space)
  dat <- dat[ind != 0.0]


  # take cumulative sum
  dat[, (dep) := cumsum((sign(get(dep)) * sign(get(ind))) * abs(get(dep)))]
  dat[, (ind) := cumsum(abs(get(ind)))]


  fit_least_squares(dat, dep, ind, return_model)

}
