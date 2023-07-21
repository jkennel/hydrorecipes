#' lag_difference
#'
#' Calculate differences for the dependent and independent variables and remove any NA values.  This modifies the data.table in place.
#'
#' @inheritParams be_least_squares_diff
#' @param var_name name of the column to lag (character)
#' @param remove_na remove NA values (logical)
#'
#' @return data.table with lagged differences
#' @importFrom stats na.omit
#' @export
#'
#' @examples
#' library(data.table)
#' datetime <- seq.POSIXt(as.POSIXct("2016-01-01 12:00:00"),
#'                        as.POSIXct("2016-01-05 12:00:00"), by='hour' )
#' baro <- sin(seq(0, 2*pi, length.out = length(datetime)))
#' wl <- 0.4 * baro
#' dat <- data.table(baro, wl, datetime)
#' lag_difference(dat, lag_space = 1, remove_na = FALSE)
#'
lag_difference <- function(dat,
                           var_name = 'wl',
                           lag_space = 1,
                           inverse = FALSE,
                           remove_na = TRUE){
  
  dat[, (var_name) := diff_shift(get(var_name), lag_space) ]
  
  
  if(remove_na){
    dat <- na.omit(dat)
  }
  
  if (inverse) {
    dat[, (var_name) := -get(var_name) ]
  }
  
  return(dat)
}


#' diff_shift
#'
#' Calculate lagged differences padded with NA values
#'
#' @param x vector to difference (numeric)
#' @param lag_space spacing for lags, useful for higher frequency monitoring
#' (integer).
#'
#' @return lagged differences
#' @export
#'
#' @examples
#' diff_shift(1:100, 2)
#'
diff_shift <- function(x,
                       lag_space=1) {
  
  return(x - data.table::shift(x, lag_space))
  
}

