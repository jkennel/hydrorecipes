#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Papadopulos-Cooper Wellbore Storage Test -------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepAquiferWellboreStorage <- R6Class(
  classname = "step_aquifer_wellbore_storage",
  inherit = Step,
  public = list(

    # step specific variables
    time = NULL,
    flow_rate = NULL,
    radius = NULL,
    radius_casing = NULL,
    radius_well = NULL,
    thickness = NULL,
    specific_storage = NULL,
    hydraulic_conductivity = NULL,
    n_terms = NULL,

    initialize = function(time,
                          flow_rate = 10.0,
                          radius = 0.15,
                          radius_casing = 0.15,
                          radius_well = 0.15,
                          thickness = 1.0,
                          specific_storage = 1.0e-6,
                          hydraulic_conductivity = 10.0,
                          n_terms = 12L,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_aquifer_wellbore_storage"
      env_list$type <- "add"
      super$initialize(
        terms = c(as.symbol(time)),
        env_list,
        ...
      )


      # step specific values
      self$time <- time
      self$flow_rate <- flow_rate
      self$radius <- radius
      self$radius_casing <- radius_casing
      self$radius_well <- radius_well
      self$thickness <- thickness
      self$specific_storage <- specific_storage
      self$hydraulic_conductivity <- hydraulic_conductivity
      self$n_terms <- n_terms

      self$columns <- time

      invisible(self)
    },
    bake = function(s) {

      self$new_columns <- self$prefix

      Tr <- self$hydraulic_conductivity * self$thickness
      S  <- self$specific_storage * self$thickness


      pc <- list(papadopulos_cooper_laplace(s[["result"]][[self$columns[[1L]]]],
                                            self$flow_rate,
                                            self$radius,
                                            self$radius_casing,
                                            self$radius_well,
                                            Tr,
                                            S,
                                            1e-8, # currently not used
                                            self$n_terms))

      self$result <- setNames(pc, self$new_columns)
      return(NULL)

    }
  )
)
