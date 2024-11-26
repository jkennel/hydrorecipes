#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Remove the Central Value (mean) from a Regressor Step ------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepBaroFrequencyUnconfined <- R6Class(
  classname = "step_baro_frequency_unconfined",
  inherit = Step,
  public = list(
    column_values = c(),
    radius_well = NULL,
    storage_aquifer = NULL,
    specific_yield = NULL,
    k_vertical = NULL,
    diffusivity_vertical = NULL,
    diffusivity_vadose = NULL,
    thickness_saturated_well = NULL,
    thickness_vadose = NULL,
    thickness_aquifer = NULL,
    loading_efficiency = NULL,
    attenuation = NULL,

    # step specific variables
    initialize = function(frequency,
                          radius_well,
                          storage_aquifer,
                          specific_yield,
                          k_vertical,
                          diffusivity_vertical,
                          diffusivity_vadose,
                          thickness_saturated_well,
                          thickness_vadose,
                          thickness_aquifer,
                          loading_efficiency,
                          attenuation,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      frequency <- deparse(substitute(frequency))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_baro_frequency_unconfined"
      env_list$type <- "add"
      super$initialize(
        terms = as.symbol(frequency),
        env_list,
        ...
      )



      self$radius_well <- radius_well
      self$storage_aquifer <- storage_aquifer
      self$specific_yield <- specific_yield
      self$k_vertical <- k_vertical
      self$diffusivity_vertical <- diffusivity_vertical
      self$diffusivity_vadose <- diffusivity_vadose
      self$thickness_saturated_well <- thickness_saturated_well
      self$thickness_vadose <- thickness_vadose
      self$thickness_aquifer <- thickness_aquifer
      self$loading_efficiency <- loading_efficiency
      self$attenuation <- attenuation

      invisible(self)
    },
    # subtract the central value from a column
    bake = function(new_data) {

      roj_1990 <- list(areal_rojstaczer_unconfined(
        new_data[[self$columns]],
        self$radius_well,
        self$storage_aquifer,
        self$specific_yield,
        self$k_vertical,
        self$diffusivity_vertical,
        self$diffusivity_vadose,
        self$thickness_saturated_well,
        self$thickness_vadose,
        self$thickness_aquifer,
        self$loading_efficiency,
        self$attenuation
      ))
      self$new_columns <- self$prefix

      self$result <- setNames(roj_1990, self$new_columns)
      self$result

    }
  )
)



