#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate Barometric Efficiency using Harmonic tide methods ------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepBaroHarmonic`
#'
#' @inheritParams Step
#' @inheritParams recipes::step_center
#'
#'
#' @references
#' Rau, G.C., Cuthbert, M.O., Acworth, R.I. and Blum, P., 2020.
#'   Disentangling the groundwater response to Earth and atmospheric tides
#'   to improve subsurface characterisation. Hydrology and earth system
#'   sciences, 24(12), pp.6033-6046.
#'
#' Acworth, R. I., Halloran, L. J. S., Rau, G. C., Cuthbert, M. O.,
#'   & Bernardi, T. L. (2016). An objective frequency-domain method for
#'   quantifying confined aquifer compressible storage using Earth and
#'   atmospheric tides. Geophysical Research Letters, 43(November).
#'   https://doi.org/10.1002/2016GL071328
#'
#' @family barometric
#'
#' @export
StepBaroHarmonic <- R6Class(
  classname = "step_baro_rau",
  inherit = Step,
  public = list(

    water_level = NULL,
    barometric_pressure = NULL,
    earth_tide = NULL,
    time = NULL,

    frequency = NULL,
    cycle_size = NULL,
    start = NULL,
    inverse = NULL,

    barometric_efficiency = list(),

    # step specific variables
    initialize = function(time,
                          water_level,
                          barometric_pressure,
                          earth_tide,
                          frequency = c(1.9324, 2.0),
                          cycle_size = 86400,
                          start = 0.0,
                          inverse = TRUE,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      water_level <- deparse(substitute(water_level))
      barometric_pressure <- deparse(substitute(barometric_pressure))
      earth_tide <- deparse(substitute(earth_tide))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_baro_harmonic"
      env_list$type <- "augment"
      super$initialize(
        terms = c(
          as.symbol(time),
          as.symbol(water_level),
          as.symbol(barometric_pressure),
          as.symbol(earth_tide)),
        env_list,
        ...
      )

      self$time <- time
      self$water_level <- water_level
      self$barometric_pressure <- barometric_pressure
      self$earth_tide <- earth_tide
      self$frequency <- frequency
      self$cycle_size <- cycle_size
      self$start <- start
      self$inverse <- inverse

      self$columns <- c(time, water_level, barometric_pressure, earth_tide)

      invisible(self)
    },
    bake = function(new_data) {

      nms <- names(new_data)

      # create regression formula
      formula_txt <- paste0(paste(nms[-1], collapse = '+'), "~", nms[1])

      # include linear trend and intercept
      harmonics <- frecipes::Recipe$new(formula = as.formula(formula_txt), new_data)$
        add_step(StepIntercept$new())$
        add_step(StepHarmonic$new(datetime,
                                  frequency = self$frequency,
                                  cycle_size = self$cycle_size,
                                  starting_value = self$start))$
        plate("m")


      X <- harmonics[, -(2:4)]
      Y <- harmonics[, (2:4)]


      soln <- llt_solve(X, Y)

      co_names <- colnames(X)
      wh_s <- grep("sin", co_names)
      wh_c <- grep("cos", co_names)
      soln_cplx <- sin_cos_to_complex(c = soln[wh_c,], s = -soln[wh_s,])

      self$barometric_efficiency <- be_harmonic_cpp(soln_cplx, self$inverse)

      print(self$barometric_efficiency)

    return(NULL)

    }
  )
)

sin_cos_to_complex <- function(c, s) {
  complex(real = t(c), imaginary = t(s))
}
