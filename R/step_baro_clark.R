#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate Barometric Efficiency using Clarks Method --------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepBaroClark <- R6Class(
  classname = "step_baro_clark",
  inherit = Step,
  public = list(

    water_level = NULL,
    barometric_pressure = NULL,
    lag_space = NULL,
    inverse = NULL,
    differences = TRUE,

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
    bake = function(s) {

      be <- list()
      for (i in seq_along(self$lag_space)) {
        be[[i]] <- list(be = be_clark_cpp(
          dep = s[["result"]][[self$columns[[1L]]]],
          ind = s[["result"]][[self$columns[[2L]]]],
          lag_space = self$lag_space[i],
          inverse = self$inverse
        ))

      }

      be <- collapse::rowbind(be)
      self$barometric_efficiency <- be

      return(NULL)
    }
  )
)
