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
    return(0:max_time_lag)
  } else {

    lags <- round(10^seq(0,
                         log10(max_time_lag + 1),
                         length.out = n) - 1)

    wh <- which(lags <= 0:(n-1))

    if(length(wh) > 0) {
      lags[wh] <- wh-1
    }

    return(as.integer(lags))
  }

}
