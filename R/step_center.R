#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Remove the Central Value (mean) from a Regressor Step ------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepCenter <- R6Class(
  classname = "step_center",
  inherit = Step,
  public = list(
    column_values = c(),
    na_rm = NA,
    fun = NULL,

    # step specific variables
    initialize = function(terms,
                          na_rm = TRUE,
                          fun = collapse::fmean,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_center"
      env_list$type <- "modify"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$na_rm <- na_rm
      self$fun <- fun

      invisible(self)
    },
    prep = function(data) {

      self$column_values <- self$fun(data, na.rm = self$na_rm)

    },

    # subtract the central value from a column
    bake = function(s) {

      s[["result"]][self$columns] <- s[["result"]][self$columns] %r-% self$column_values

      return(NULL)

    }
  )
)
