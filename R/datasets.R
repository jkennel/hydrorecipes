#' @title transducer
#'
#' @description This `data.frame` contains the water levels, barometric pressure,
#' and synthetic Earth tides from 2016-08-25 to 2016-10-15 (sub-sampled to every
#' 2 minutes).
#'
#' @format A \code{data.frame} The columns are:
#' \describe{
#'  \item{\code{datetime}}{POSIXct date and time}
#'  \item{\code{wl}}{transducer pressure of water column and air (dbar)}
#'  \item{\code{baro}}{barometric pressure (dbar)}
#'  \item{\code{et}}{synthetic gravity}
#' }
#'
#' @examples
#' data(transducer)
'transducer'


#' @title wipp30
#'
#' @description This `data.frame` contains the water levels, barometric pressure,
#' and Earth tides for wipp30.  This is the dataset included in BETCO.
#'
#' @format A \code{data.frame} The columns are:
#' \describe{
#'  \item{\code{time}}{elapsed time in hours}
#'  \item{\code{wl}}{water level}
#'  \item{\code{baro}}{barometric pressure}
#'  \item{\code{et}}{synthetic gravity}
#' }
#'
#' @references
#' Toll NJ, Rasmussen TC, 2007, "Removal of barometric pressure effects and
#'   Earth tides from observed water levels", Ground Water 45(1):101-105
#'
#' @examples
#' data(wipp30)
'wipp30'
