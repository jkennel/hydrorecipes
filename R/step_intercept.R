#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Add an Intercept Term --------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepIntercept <- R6Class(
  classname = "step_intercept",
  inherit = Step,
  public = list(

    value = NA_real_,
    # step specific variables
    initialize = function(terms,
                          value = 1.0,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      env_list <- get_function_arguments()
      env_list$step_name <- "step_intercept"
      env_list$type <- "add"
      super$initialize(
        terms = NULL,
        env_list[names(env_list) != "terms"],
        ...
      )
      self$value <- value


      invisible(self)
    },
    bake = function(s) {
      n <- length(s[["result"]][[1L]])
      self$new_columns <- self$prefix

      self$result <- setNames(list(rep.int(self$value, n)), self$new_columns)
      return(NULL)


    }
  )
)
