#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Remove the Central Value (mean) from a Regressor Step ------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepTidalCooper <- R6Class(
  classname = "step_tidal_cooper",
  inherit = Step,
  public = list(
    column_values = c(),
    frequency = NULL,
    storativity = NULL,
    transmissivity = NULL,
    thickness_aquifer = NULL,
    height_water = NULL,
    radius_well = NULL,
    radius_casing = NULL,
    gravity =  9.80665,

    # step specific variables
    initialize = function(frequency,
                          storativity,
                          transmissivity,
                          thickness_aquifer,
                          height_water,
                          radius_well,
                          radius_casing = radius_well,
                          gravity =  9.80665,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      frequency <- deparse(substitute(frequency))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_tidal_cooper"
      env_list$type <- "add"
      super$initialize(
        terms = as.symbol(frequency),
        env_list,
        ...
      )

      self$storativity = storativity
      self$transmissivity = transmissivity
      self$thickness_aquifer = thickness_aquifer
      self$height_water = height_water
      self$radius_well = radius_well
      self$radius_casing = radius_casing
      self$gravity =  gravity

      invisible(self)
    },
    # subtract the central value from a column
    bake = function(new_data) {

      cooper_1965 <- list(
        tidal_cooper_1965(
          self$frequency,
          self$storativity,
          self$transmissivity,
          self$thickness_aquifer,
          self$height_water,
          self$radius_well,
          self$radius_casing,
          self$gravity))

      self$new_columns <- self$prefix

      setNames(cooper_1965, self$new_columns)

    }
  )
)
