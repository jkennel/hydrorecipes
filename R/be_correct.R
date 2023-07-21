#' be_correct
#'
#' Adjust values based on the barometric efficiency
#'
#' @inheritParams be_least_squares_diff
#' @param be \code{numeric} value of the barometric efficiency
#' @param known_mean \code{numeric} explicitly enter the mean if known.  Otherwise estimate from data.
#'
#' @return \code{numeric vector} of corrected values
#' @export
#'
#' @examples
#' library(data.table)
#' baro = rnorm(1000, sd = 0.01) + 9
#' wl <- baro * 0.4 + 18
#' dat <- data.table(baro, wl)
#' dat$wl + be_correct(dat, be=0.4, inverse = FALSE)
#' dat$wl + be_correct(dat, be=0.4, inverse = FALSE, known_mean = 9) 
#' # should return ~21.6 = 18 + 0.4 * 9
#'
be_correct <- function(dat,
                       dep = 'wl',
                       ind = 'baro',
                       be = 0,
                       inverse = TRUE,
                       known_mean = NULL){
  
  # hack for 'global variables NOTE
  corrected <- NULL
  
  dat_c <- copy(dat)
  adj   <- dat_c[[ind]]
  
  
  adj[is.na(adj)] <- 0
  if (is.null(known_mean)){
    adj <- be * (adj - mean(adj, na.rm = TRUE))
  } else {
    adj <- be * (adj - known_mean)
  }
  
  if (inverse) {
    return(adj)
  } else {
    return(-adj)
  }
  
}


