#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Remove Regressors Step -------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepDropColumns <- R6Class(
  classname = "step_drop_columns",
  inherit = Step,
  public = list(
    initialize = function(terms,
                          role = "modify",
                          ...) {

      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_drop_columns"
      env_list$type <- "modify"
      super$initialize(terms = terms,
                       env_list[names(env_list) != "terms"],
                       ...
      )


      invisible(self)
    },

    bake = function(s) {

      s[["result"]][self$columns] <- NULL

      return(NULL)
    }

  )
)
