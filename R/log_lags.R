#' @title
#' log_lags
#'
#' @description
#' Generate lags or knots with logarithmic spacing. Lags/knots start at 0. Lags
#' are in terms of samples. For example, if samples are taken every 60, a
#' lag of 1 corresponds to 60 seconds, a lag of 4 would correspond to 240
#' seconds.
#'
#' @param n The total number of lags (integer)
#' @param max_time_lag The maximum lag in number of samples (integer)
#'
#' @return An integer vector of lags
#'
#' @export
#'
#' @examples
#' log_lags(12, 86400)
log_lags <- function(n, max_time_lag) {

  if(length(n) != 1) {
    stop('n must be length 1')
  }

  if(length(max_time_lag) != 1) {
    stop('n must be length 1')
  }

  if(n <= 0) {
    stop('n must be greater than 0')
  }

  if(max_time_lag < 0) {
    stop('max_time_lag must be non-negative')
  }

  if(!is.numeric(n)) {
    stop('n must be an integer or coercible to an integer')
  }

  if(!is.numeric(max_time_lag)) {
    stop('max_time_lag must be an integer or coercible to an integer')
  }

  if(!is.integer(n)) {
    n <- as.integer(n)
  }

  if(!is.integer(max_time_lag)) {
    max_time_lag <- as.integer(max_time_lag)
  }

  if(n > (max_time_lag + 1L)) {
    warning('The number of lags is greater than the maximum time lag')
    return(0L:max_time_lag)
  }

  lags <- round(10^seq(0L,
                       log10(max_time_lag + 1L),
                       length.out = n) - 1L)

  wh <- which(lags <= 0L:(n-1L))

  if(length(wh) > 0L) {
    lags[wh] <- wh - 1L
  }

  as.integer(lags)

}
