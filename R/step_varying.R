#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Remove Columns with a Single Value -------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepVarying <- R6Class(
  classname = "step_varying",
  inherit = Step,
  public = list(
    to_remove = NULL,

    # step specific variables
    initialize = function(terms,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_varying"
      env_list$type <- "modify"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )


      invisible(self)
    },
    bake = function(new_data) {
      self$to_remove <- !collapse::varying(new_data)
      new_data[self$to_remove] <- list(NULL)

      self$result <- new_data
      self$result
    }
  )
)
