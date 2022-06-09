#' @title
#' log_lags
#'
#' @description
#' Generate lags that are spaced logarithmic. Lags start at 0.
#'
#' @param n number of lags (integer)
#' @param max_time_lag maximum lag (integer)
#'
#' @return integer vector of lags
#'
#' @export
#'
#' @examples
#' log_lags(12, 86401)
log_lags <- function(n, max_time_lag) {


  if(n <= 0) {
    stop('n must be greater than 0')
  }


  if(!is.numeric(n)) {
    stop('n must be an integer or coercible to an integer')
  }

  if(!is.numeric(max_time_lag)) {
    stop('max_time_lag must be an integer or coercible to an integer')
  }


  if(!is.integer(n) | !is.integer(max_time_lag)) {
    n <- as.integer(n)
    max_time_lag <- as.integer(max_time_lag)
  }


  if(n > max_time_lag) {
    warning('The number of lags is greater than the maximum time lag')
    return(0L:max_time_lag)
  }

  lags <- round(10^seq(0L,
                       log10(max_time_lag + 1L),
                       length.out = n) - 1L)

  wh <- which(lags <= 0L:(n-1L))

  if(length(wh) > 0L) {
    lags[wh] <- wh-1
  }

  as.integer(lags)

}
