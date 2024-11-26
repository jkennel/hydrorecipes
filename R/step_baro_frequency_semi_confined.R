#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Remove the Central Value (mean) from a Regressor Step ------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepBaroFrequencySemiConfined <- R6Class(
  classname = "step_baro_frequency_semi_confined",
  inherit = Step,
  public = list(
    column_values = c(),
    radius_well = NULL,
    transmissivity = NULL,
    storage_confining = NULL,
    storage_aquifer = NULL,
    diffusivity_confining = NULL,
    diffusivity_vadose = NULL,
    thickness_confining = NULL,
    thickness_vadose = NULL,
    loading_efficiency = NULL,
    attenuation = NULL,

    # step specific variables
    initialize = function(frequency,
                          radius_well,
                          transmissivity,
                          storage_confining,
                          storage_aquifer,
                          diffusivity_confining,
                          diffusivity_vadose,
                          thickness_confining,
                          thickness_vadose,
                          loading_efficiency,
                          attenuation,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      frequency <- deparse(substitute(frequency))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_baro_frequency_semi_confined"
      env_list$type <- "add"
      super$initialize(
        terms = as.symbol(frequency),
        env_list,
        ...
      )



      self$radius_well <- radius_well
      self$transmissivity <- transmissivity
      self$storage_confining <- storage_confining
      self$storage_aquifer <- storage_aquifer
      self$diffusivity_confining <- diffusivity_confining
      self$diffusivity_vadose <- diffusivity_vadose
      self$thickness_confining <- thickness_confining
      self$thickness_vadose <- thickness_vadose
      self$loading_efficiency <- loading_efficiency
      self$attenuation <- attenuation

      invisible(self)
    },
    # subtract the central value from a column
    bake = function(new_data) {

      roj_1988 <- list(areal_rojstaczer_semiconfined(
        new_data[[self$columns]],
        self$radius_well,
        self$transmissivity,
        self$storage_confining,
        self$storage_aquifer,
        self$diffusivity_confining,
        self$diffusivity_vadose,
        self$thickness_confining,
        self$thickness_vadose,
        self$loading_efficiency,
        self$attenuation
      ))

      self$new_columns <- self$prefix

      self$result <- setNames(roj_1988, self$new_columns)
      self$result

    }
  )
)
