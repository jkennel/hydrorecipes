#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Theis Step -------------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepAquiferTheisAniso <- R6Class(
  classname = "step_aquifer_theis_aniso",
  inherit = Step,



  public = list(
    # step specific variables
    time = NULL,
    flow_rate = NULL,
    thickness = NULL,
    distance_x = NULL,
    distance_y = NULL,
    specific_storage = NULL,
    hydraulic_conductivity_x = NULL,
    hydraulic_conductivity_y = NULL,

    initialize = function(time,
                          flow_rate,
                          thickness = 1.0,
                          distance_x = 100.0,
                          distance_y = 100.0,
                          specific_storage = 1.0e-6,
                          hydraulic_conductivity_x = 1.0e-4,
                          hydraulic_conductivity_y = 1.0e-4,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      flow_rate <- deparse(substitute(water_level))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_aquifer_theis_aniso"
      env_list$type <- "predictor"
      super$initialize(
        terms = c(
          as.symbol(time),
          as.symbol(flow_rate)),
        env_list,
        ...
      )

      self$time <- time
      self$flow_rate <- flow_rate
      self$distance_x <- distance_x
      self$distance_y <- distance_y
      self$specific_storage <- specific_storage
      self$hydraulic_conductivity_x <- hydraulic_conductivity_x
      self$hydraulic_conductivity_y <- hydraulic_conductivity_y

      self$columns <- c(time, flow_rate)

      invisible(self)
    },
    bake = function(new_data) {

      self$new_columns <- self$prefix
      self$columns <- paste(self$columns, collapse = ",")

      setNames(theis_aniso_time(
        distance_x = self$distance_x,
        distance_y = self$distance_y,
        specific_storage = self$specific_storage,
        hydraulic_conductivity_x = self$hydraulic_conductivity,
        hydraulic_conductivity_y = self$hydraulic_conductivity_y,
        thickness = self$thickness,
        time = new_data[[self$time]],
        flow_rate = new_data[[self$flow_rate]]
      ), self$new_columns)


    }
  )
)


