#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Vary Step -------------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepVaryParameter <- R6Class(
  classname = "step_vary_parameter",
  inherit = Step,

  public = list(

    # step specific variables
    initialize = function(parameters,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      env_list <- get_function_arguments()
      env_list$step_name <- "step_vary_parameter"
      env_list$type <- "modify"
      super$initialize(
        terms = NULL,
        env_list[names(env_list) != "terms"],
        ...
      )

      do.call(super$initialize, inputs)
      self$parameters <- parameters

      invisible(self)

    }
  )
)
