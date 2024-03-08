#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate Barometric Efficiency using Clarks Method --------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepBaroClark` adjust the central value to zero.
#'
#' @inheritParams Step
#' @inheritParams recipes::step_center
#'
#'
#' @family barometric
#'
#' @export
StepBaroClark <- R6Class(
  classname = "step_baro_clark",
  inherit = Step,
  public = list(

    water_level = NULL,
    barometric_pressure = NULL,
    lag_space = NULL,
    inverse = NULL,

    barometric_efficiency = c(),

    # step specific variables
    initialize = function(water_level,
                          barometric_pressure,
                          lag_space = 1L,
                          inverse = FALSE,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      water_level <- deparse(substitute(water_level))
      barometric_pressure <- deparse(substitute(barometric_pressure))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_baro_clark"
      env_list$type <- "augment"
      super$initialize(
        terms = c(as.symbol(water_level), as.symbol(barometric_pressure)),
        env_list,
        ...
      )

      self$water_level <- water_level
      self$barometric_pressure <- barometric_pressure
      self$lag_space <- lag_space
      self$inverse <- inverse
      self$columns <- c(water_level, barometric_pressure)

      invisible(self)
    },
    bake = function(new_data) {

      for (i in seq_along(self$lag_space)) {
        self$barometric_efficiency[i] <- be_clark_cpp(
          dep = new_data[[1]],
          ind = new_data[[2]],
          lag_space = self$lag_space[i],
          inverse = self$inverse
        )
      }

    return(NULL)
    }
  )
)
