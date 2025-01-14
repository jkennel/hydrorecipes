#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate Barometric Efficiency using Least Squares --------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepBaroLeastSquares <- R6Class(
  classname = "step_baro_least_squares",
  inherit = Step,
  public = list(

    water_level = NULL,
    barometric_pressure = NULL,
    lag_space = NULL,
    inverse = NULL,
    differences = NULL,

    barometric_efficiency = list(),

    # step specific variables
    initialize = function(water_level,
                          barometric_pressure,
                          lag_space = 1L,
                          inverse = FALSE,
                          differences = FALSE,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      water_level <- deparse(substitute(water_level))
      barometric_pressure <- deparse(substitute(barometric_pressure))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_baro_least_squares"
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
      self$differences <- differences
      self$columns <- c(water_level, barometric_pressure)

      invisible(self)
    },
    bake = function(s) {

      be <- list()
      for (i in seq_along(self$lag_space)) {
        if (self$differences) {
          be[[i]] <- list(be = be_least_squares_diff_cpp(
            dep = s[["result"]][[self$columns[[1L]]]],
            ind = s[["result"]][[self$columns[[2L]]]],
            lag_space = self$lag_space[i],
            inverse = self$inverse
          ))
        } else {

          be[[i]] <- list(be = be_least_squares_cpp(
            dep = s[["result"]][[self$columns[[1L]]]],
            ind = s[["result"]][[self$columns[[2L]]]],
            inverse = self$inverse
          ))

        }
        be <- collapse::rowbind(be)
        self$barometric_efficiency <- be
      }

      return(NULL)
    }
  )
)
