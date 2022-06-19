#' @title transducer
#'
#' @description This data.frame contains the water levels, barometric pressure,
#' and synthetic Earth tides from 2016-08-25 to 2016-10-15 (sub-sampled to every
#' 2 minutes).
#'
#' @format A \code{data.frame} The columns are:
#' \describe{
#'  \item{\code{datetime}}{POSIXct date and time}
#'  \item{\code{wl}}{transducer pressure of water column(dbar)}
#'  \item{\code{baro}}{barometric pressure (dbar)}
#'  \item{\code{et}}{synthetic gravity}
#' }
#'
#' @examples
#' utils::data(transducer)
'transducer'


#' @title wipp30
#'
#' @description This data.frame contains the water levels, barometric pressure,
#' and Earth tides for wipp30.  This is the dataset included in BETCO.
#'
#' @format A \code{data.frame} The columns are:
#' \describe{
#'  \item{\code{time}}{time in hours}
#'  \item{\code{wl}}{transducer pressure of water column}
#'  \item{\code{baro}}{barometric pressure}
#'  \item{\code{et}}{synthetic gravity}
#' }
#'
#' @examples
#' utils::data(wipp30)
'wipp30'
