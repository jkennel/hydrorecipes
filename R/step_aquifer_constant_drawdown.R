#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Jacob-Lohman Constant Drawdown Test ------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepAquiferConstantDrawdown <- R6Class(
  classname = "step_aquifer_constant_drawdown",
  inherit = Step,
  public = list(

    # step specific variables
    time = NULL,
    drawdown = NULL,
    thickness = NULL,
    radius_well = NULL,
    specific_storage = NULL,
    hydraulic_conductivity = NULL,
    n_terms = NULL,

    initialize = function(time,
                          drawdown = 1.0,
                          thickness = 1.0,
                          radius_well = 0.15,
                          specific_storage = 1.0e-6,
                          hydraulic_conductivity = 1.0e-4,
                          n_terms = 16L,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_aquifer_constant_drawdown"
      env_list$type <- "add"
      super$initialize(
        terms = as.symbol(time),
        env_list,
        ...
      )


      # step specific values
      self$time <- time
      self$drawdown <- drawdown
      self$thickness <- thickness
      self$radius_well <- radius_well
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


      # add result data
      self$result <- setNames(list(jacob_lohman_laplace(
        time = s[["result"]][[self$columns]],
        s = self$drawdown,
        r = self$radius_well,
        Tr = Tr,
        S = S,
        prec = 1e-8, # not currently used
        n_terms = self$n_terms
      )), self$new_columns)


      return(NULL)
    }
  )
)
