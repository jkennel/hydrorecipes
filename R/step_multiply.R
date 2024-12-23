#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Adjust the dispersion (e.g. scale by standard deviation) ---------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepMultiply <- R6Class(
  classname = "step_multiply",
  inherit = Step,
  public = list(
    values = NA_real_,

    # step specific variables
    initialize = function(terms,
                          values = 1.0,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_multiply"
      env_list$type <- "modify"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$values <- values

      invisible(self)
    },

    # multiply column
    bake = function(s) {

      self$result <- s[["result"]][self$columns] %r*% self$values

      return(NULL)

    }

  )
)
