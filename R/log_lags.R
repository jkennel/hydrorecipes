#' @title
#' log_lags
#'
#' @description
#' Generate logarithmically spacing for lags. Note: the number of lags will not
#' exactly equal n unless max_time_lag is large or n is very small
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

  if(n > max_time_lag) {
    warning('The number of lags is greater than the maximum time lag')
    return(1:max_time_lag - 1)
  } else {
    lags <- round(10^seq(0,
                   log10(max_time_lag),
                   length.out = n) - 1)

    wh <- which(lags <= 1:n)

    if(length(wh) > 0) {
      lags[wh] <- wh - 1
    }

    return(lags)
  }

}
