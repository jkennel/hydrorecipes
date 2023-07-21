#' be_ratio
#'
#' Calculate the barometric efficiency using the mean (or other statistic) 
#' of the water level change divided by the barometric pressure change .  
#' There is the option to only include responses that are large.
#'
#' @inheritParams be_least_squares_diff
#' @param quant quantile cutoff value that differences need to be this large (numeric)
#' @param stat the result to return (mean, median, quantile) (function)
#' @param ... other arguments to pass to "stat"
#'
#' @return barometric efficiency calculated by using ratio method
#' @export
#'
#' @examples
#' library(data.table)
#' datetime <- seq.POSIXt(as.POSIXct("2016-01-01 12:00:00"),
#'                        as.POSIXct("2016-01-05 12:00:00"), by='hour' )
#' baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
#' noise <- rnorm(length(datetime), sd = 0.01)
#' wl <- -0.4 * baro + noise
#' dat <- data.table(baro, wl, datetime)
#' be_ratio(dat, quant = 0.5, stat=median)
#' be_ratio(dat, quant = 0.9, stat=mean)
#'
be_ratio <- function(dat,
                     dep = 'wl',
                     ind = 'baro',
                     lag_space = 1,
                     inverse = TRUE,
                     quant = 0.90,
                     stat = mean,
                     ...){
  
  # hack for 'global variables NOTE
  ratio <- NULL
  
  dat <- dat[, c(dep, ind), with = FALSE]
  
  dat <- lag_difference(dat, dep, lag_space, inverse)
  dat <- lag_difference(dat, ind, lag_space)
  
  dat <- dat[get(ind) != 0]
  
  # select the larger changes in barometric pressure
  dat <- dat[abs(get(ind)) >= stats::quantile(abs(get(ind)), quant, na.rm = TRUE) ]
  
  dat[, ratio := get(dep) / get(ind) ]
  
  # minimize the effect of wl changes that are not a result of barometric changes
  # be values should be between 0.0 and 1.0
  # dat <- dat[ ratio >= 1, ratio:=1 ]
  # dat <- dat[ ratio <= 0, ratio:=0 ]
  
  stat(dat$ratio, na.rm = TRUE, ...)
  
}

