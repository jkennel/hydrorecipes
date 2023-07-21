#' @title
#' be_rahi ## need to check for correctness
#'
#' @description
#' Rahi 2010 solution for calculating barometric efficiency.
#'
#' @inheritParams be_least_squares_diff
#'
#' @return barometric efficiency using Rahi's method
#' 
#' @references Rahi, K. A. (2010). Estimating the hydraulic parameters of the 
#' Arbuckle-Simpson aquifer by analysis of naturally-induced stresses 
#' (Doctoral dissertation, Oklahoma State University).
#' 
#' @export
#'
#' @examples
#' library(data.table)
#' datetime <- seq.POSIXt(as.POSIXct("2016-01-01 12:00:00"),
#'                        as.POSIXct("2016-01-05 12:00:00"), by='hour' )
#' baro <- sin(seq(0, 2 * pi, length.out = length(datetime)))
#' noise <- rnorm(length(datetime), sd = 0.01)
#' wl <- -0.4 * baro + noise
#' dat <- data.table(baro, wl, datetime)
#' be_rahi(dat, dep = 'wl', ind = 'baro', lag_space = 1)
#'
be_rahi <- function(dat,
                    dep = 'wl',
                    ind = 'baro',
                    lag_space = 1,
                    inverse = TRUE) {
  
  
  # subset the necessary variables, this also does a copy
  dat <- dat[, c(dep, ind), with = FALSE]
  
  
  # do first difference
  dat <- lag_difference(dat, dep, lag_space, inverse)
  dat <- lag_difference(dat, ind, lag_space)
  
  # don't count if the signs are different
  dat <- dat[sign(get(dep)) != sign(get(ind)), c(dep, ind) := 0.0]
  
  # water level change must be less than baro change
  dat[abs(get(dep)) > abs(get(ind)), c(dep, ind) := 0.0]
  
  
  # return the ratio
  return(sum(abs(dat[[dep]])) / sum(abs(dat[[ind]])))
  
}
