#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate Barometric Efficiency using Acworth Method -------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepBaroAcworth`
#'
#' @inheritParams Step
#' @inheritParams recipes::step_center
#'
#'
#' @family barometric
#'
#' @export
StepBaroAcworth <- R6Class(
  classname = "step_baro_acworth",
  inherit = Step,
  public = list(

    water_level = NULL,
    barometric_pressure = NULL,
    earth_tides = NULL,

    frequency_a = NULL,
    frequency_b = NULL,
    inverse = NULL,
    spans = NULL,
    detrend = NULL,
    demean = NULL,
    taper = NULL,


    barometric_efficiency = c(),

    # step specific variables
    initialize = function(water_level,
                          barometric_pressure,
                          earth_tides,
                          frequency_a = 1.9324, # m2
                          frequency_b = 2.0,    # s2
                          inverse = FALSE,
                          spans = 5,
                          detrend = TRUE,
                          demean = TRUE,
                          taper = 0.1,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      water_level <- deparse(substitute(water_level))
      barometric_pressure <- deparse(substitute(barometric_pressure))
      earth_tides <- deparse(substitute(earth_tides))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_baro_acworth"
      env_list$type <- "augment"
      super$initialize(
        terms = c(as.symbol(water_level),
                  as.symbol(barometric_pressure),
                  as.symbol(earth_tides)),
        env_list,
        ...
      )

      self$water_level <- water_level
      self$barometric_pressure <- barometric_pressure
      self$earth_tides <- earth_tides
      self$frequency_a <- frequency_a
      self$frequency_b <- frequency_b
      self$inverse <- inverse
      self$spans <- spans
      self$detrend <- detrend
      self$demean <- demean
      self$taper <- taper
      self$columns <- c(water_level, barometric_pressure, earth_tides)

      invisible(self)
    },
    bake = function(new_data) {

      self$barometric_efficiency <- be_acworth_cpp(
        x = collapse::qM(new_data),
        spans = self$spans,
        detrend = self$detrend,
        demean = self$demean,
        taper = self$taper,
        inverse = self$inverse,
        f1 = self$frequency_a,
        f2 = self$frequency_b,
        frequency_scale = self$frequency_scale
      )

    return(NULL)
    }
  )
)
