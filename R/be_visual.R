#' be_visual_data
#'
#' Generate dataset for comparing barometric efficiency
#'
#' @inheritParams be_visual
#'
#' @return data.table of barometric efficiency compensated datasets
#' @export
#'
#' @examples
#' library(data.table)
#' be <- 0.43
#' x <- seq(0, 28*pi, pi / (12*12))
#'
#' baro <- sin(x) + rnorm(length(x), sd = 0.04)
#' wl <- -sin(x) * be + rnorm(length(x), sd = 0.04)
#' dat <- data.table(datetime = as.POSIXct(x * 86400 / (2 * pi),
#'                                         origin = '1970-01-01', tz = 'UTC'),
#'                   wl = wl, baro = baro)
#' be_visual_data(dat)
#'
be_visual_data <- function(dat,
                           dep = 'wl',
                           ind = 'baro',
                           be_tests = seq(0, 1, 0.1),
                           inverse = TRUE) {

  # hack for 'global variables NOTE
  correction <- NULL
  corrected <- NULL
  be <- NULL

  dat_tmp <- copy(dat)


  dat_list <- list()
  for (i in seq_along(be_tests)) {

    dat_tmp[, correction := be_correct(dat, dep, ind, be_tests[i], inverse = inverse) ]
    dat_tmp[, corrected := get(dep) + correction]
    dat_tmp[, be := as.numeric(be_tests[i])]

    dat_list[[i]] <- copy(dat_tmp)
  }


  return(rbindlist(dat_list))

}

#' be_visual_plot
#'
#' Plot to compare barometric efficiency. Large datasets may take a long time to plot.  Subsample should be set to TRUE
#' @inheritParams be_visual
#'
#' @param subsample should the data be subsampled for plotting? (logical)
#'
#' @return plotly graph for barometric efficiency estimation with Smith method
#' @export
#'
#' @importFrom plotly plot_ly add_lines layout hide_legend animation_opts
#' @importFrom plotly animation_slider
#'
#' @examples
#' library(plotly)
#' library(data.table)
#' be <- 0.43
#' x <- seq(0, 28 * pi, pi / (12 * 12))
#'
#' baro <- sin(x) + rnorm(length(x), sd = 0.04)
#' wl <- -sin(x) * be + rnorm(length(x), sd = 0.04)
#' dat <- data.table(datetime = as.POSIXct(x * 86400 / (2 * pi),
#'                                         origin = '1970-01-01', tz = 'UTC'),
#'                   wl = wl, baro = baro)
#' dat_be <- be_visual_data(dat)
#' #be_visual_plot(dat_be) #not run
#'
be_visual_plot <- function(dat,
                           time = 'datetime',
                           subsample = TRUE){


  if (time == 'time'){
    dat[, datetime := time]
  } else if (time != 'time') {

    if (time == 'datetime') {
    } else if ('time' %in% names(dat)) {
      dat <- dat[, -c('time'), with = FALSE]
      dat[, datetime := get(time)]
    }
  }


  # hack for 'global variables NOTE
  corrected <- NULL
  datetime <- NULL
  be <- NULL


  setkey(dat, be, datetime)


  if (nrow(dat) > 60000) {

    n_group <- length(unique(dat$be))
    n <- nrow(dat)/n_group
    max_per_group <- round(60000/n_group)
    dat <- dat[, .SD[seq(1, n, length.out = max_per_group)], by = be]

  }

  # return plotly plot
  p1 <- plot_ly(dat, x = ~datetime, y = ~corrected, height = 400, width = 700)
  p1 <- add_lines(p1, frame = ~be)
  p1 <- layout(p1, xaxis = list(range = range(dat$datetime),
                                        title = ''),
                       yaxis = list(title = 'BE compensated water level'))

  p1 <- hide_legend(p1)
  p1 <- animation_opts(p1, transition = 0)
  p1 <- animation_slider(p1, currentvalue = list(prefix = "BE: ", font = list(color = "steelblue")))
  p1

}


#' be_visual
#'
#' Generate dataset for comparing barometric efficiency
#'
#' @inheritParams be_least_squares_diff
#' @param time name of the column containing the time (character)
#' @param be_tests vector of barometric efficiencies to test (between 0 and 1) (numeric)
#' @param inverse  whether the barometric relationship is inverse
#' (TRUE means that when the barometric pressure goes up the measured water
#' level goes down (vented transducer, depth to water), FALSE means that when
#' the barometric pressure goes up so does the measured pressure
#' (non-vented transducer)) (logical).
#' @param subsample should the data be subsampled for plotting? (logical)
#'
#' @return data.table of barometric efficiency compensated datasets
#'
#' @references Smith, L. A., van der Kamp, G., & Hendry, M. J. (2013). A new
#'  technique for obtaining high‐resolution pore pressure records in thick
#'  claystone aquitards and its use to determine in situ compressibility.
#'  Water Resources Research, 49(2), 732-743.
#'  \doi{https://doi.org/10.1002/wrcr.20084}
#'
#' @export
#'
#' @examples
#' library(data.table)
#' be <- 0.43
#' x <- seq(0, 28*pi, pi / (12*12))
#'
#' baro <- sin(x) + rnorm(length(x), sd = 0.04)
#' wl <- -sin(x) * be + rnorm(length(x), sd = 0.04)
#' dat <- data.table(datetime = as.POSIXct(x * 86400 / (2 * pi),
#'                                         origin = '1970-01-01', tz = 'UTC'),
#'                   wl = wl, baro = baro)
#' be_visual(dat)
#'
be_visual <- function(dat,
                      dep = 'wl',
                      ind = 'baro',
                      time = 'datetime',
                      be_tests = seq(0, 1, 0.1),
                      inverse = TRUE,
                      subsample = TRUE) {


  dt <- be_visual_data(dat, dep, ind, be_tests, inverse)
  be_visual_plot(dt, time, subsample)

}

