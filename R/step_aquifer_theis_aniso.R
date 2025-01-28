#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Theis Anisotropic Step -------------------------------------------------------
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
    storativity = NULL,
    transmissivity_major = NULL,
    transmissivity_minor = NULL,
    anisotropy = NULL,
    major_axis_angle = NULL,
    distance_x_transformed = NULL,
    distance_y_transformed = NULL,


    initialize = function(time,
                          flow_rate,
                          thickness = 1.0,
                          distance_x = 100.0,
                          distance_y = 100.0,
                          specific_storage = 1.0e-6,
                          hydraulic_conductivity_major = 1.0e-4,
                          hydraulic_conductivity_minor = 1.0e-5,
                          major_axis_angle = 0.0,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      flow_rate <- deparse(substitute(flow_rate))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_aquifer_theis_aniso"
      env_list$type <- "add"
      super$initialize(
        terms = c(
          as.symbol(time),
          as.symbol(flow_rate)),
        env_list,
        ...
      )

      # rotate coordinates to align with major axis
      cr <- coordinate_rotate(matrix(c(distance_x, distance_y), ncol = 2),
                              major_axis_angle = major_axis_angle)

      self$time <- time
      self$flow_rate <- flow_rate
      self$distance_x_transformed <- cr[1L]
      self$distance_y_transformed <- cr[2L]
      self$storativity <- specific_storage * thickness
      self$transmissivity_major <- hydraulic_conductivity_major * thickness
      self$transmissivity_minor <- hydraulic_conductivity_minor * thickness

      self$thickness <- 1.0
      self$anisotropy = hydraulic_conductivity_major / hydraulic_conductivity_minor

      self$columns <- c(time, flow_rate)

      invisible(self)
    },
    bake = function(s) {

      self$new_columns <- self$prefix
      self$columns <- paste(self$columns, collapse = ",")

      self$result <- setNames(
        theis_aniso_time(
          distance_x = self$distance_x_transformed,
          distance_y = self$distance_y_transformed,
          storativity = self$storativity,
          transmissivity_x = self$transmissivity_major,
          transmissivity_y = self$transmissivity_minor,
          thickness = self$thickness,
          time = s[["result"]][[self$columns[[1L]]]],
          flow_rate = s[["result"]][[self$columns[[2L]]]]
        ), self$new_columns)

      return(NULL)
    }

  )
)


