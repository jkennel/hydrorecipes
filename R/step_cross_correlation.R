#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Do autocorrelation or cross-correlation ------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepCrossCorrelation <- R6Class(
  classname = "step_center",
  inherit = Step,
  public = list(
    lag_max = NULL,
    acf = list(),

    # step specific variables
    initialize = function(terms,
                          lag_max = 100,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_cross_correlation"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$lag_max <- lag_max

      invisible(self)
    },

    bake = function(s) {

      if (length(self$columns) == 1) {
        self$columns[2] <- self$columns[1]
      }

      self$new_columns <- paste(self$prefix,
                                self$columns[1],
                                self$columns[2],
                                sep = "_")

      self$acf <- setNames(list(
        convolve_correlation(s[["result"]][[self$columns[1]]],
                             s[["result"]][[self$columns[2]]],
                             self$lag_max)),
        self$new_columns)

      return(NULL)

    }
  )
)
