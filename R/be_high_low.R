#' be_high_low
#'
#' Uses the maximum and minimum values for a time period to calculate be.
#'
#' @inheritParams be_least_squares_diff
#' @param time name of the column containing the time (character)
#' @param time_group format for grouping time (ie: "\%Y-\%m-\%d") (character)
#'
#' @return barometric efficiency based on high and low values
#' @export
#'
#' @examples
#' library(data.table)
#' datetime <- seq.POSIXt(as.POSIXct("2016-01-01 00:00:00", tz = 'UTC'),
#'                        as.POSIXct("2016-01-06 00:00:00", tz = 'UTC'), by='hour')
#' baro  <- sin(seq(0, 2*pi, length.out = length(datetime)))
#' wl    <- 0.4 * baro + rnorm(length(datetime), sd = 0.01)
#' dat   <- data.table(baro, wl, datetime)
#' be_hl <- be_high_low(dat)
#' be_hl <- be_high_low(dat, time_group = 86400)
#'
be_high_low <- function(dat,
                        dep = 'wl',
                        ind = 'baro',
                        time = 'datetime',
                        time_group = "%Y-%m-%d"){
  
  # hack for 'global variables NOTE
  diff_wl   <- NULL
  diff_baro <- NULL
  ratio     <- NULL
  group     <- NULL
  
  # subset the necessary variables
  dat <- na.omit(dat[, c(time, dep, ind), with = FALSE])
  
  
  if(is.numeric(time_group)){
    
    # calculate high and low values according to time grouping
    dat <- dat[, list(
      start = min(get(time)),
      end   = max(get(time)),
      diff_wl = max(get(dep)) - min(get(dep)),
      diff_baro = max(get(ind)) - min(get(ind)),
      n = .N),
      by = list(group = floor(as.numeric(get(time))/time_group))]
    
  } else {
    # calculate high and low values according to time grouping
    dat <- dat[, list(
      start = min(get(time)),
      end   = max(get(time)),
      diff_wl = max(get(dep)) - min(get(dep)),
      diff_baro = max(get(ind)) - min(get(ind)),
      n = .N),
      by = list(group = format(get(time), format = time_group))]
  }
  
  
  dat[, group := NULL]
  dat[, ratio := diff_wl / diff_baro]
  
  dat
}

