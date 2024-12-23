#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Add noise to variable --------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepAddNoise <- R6Class(
  classname = "step_add_noise",
  inherit = Step,
  public = list(

    sd = NULL,
    mean = NULL,
    fun = NULL,

    # step specific variables
    initialize = function(terms,
                          sd = 1.0,
                          mean = 0.0,
                          fun = rnorm,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_add_noise"
      env_list$type <- "modify"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$mean <- mean
      self$sd <- sd
      self$fun <- fun

      invisible(self)
    },

    # subtract the central value from a column
    bake = function(s) {

      n <- length(s[["result"]][[self$columns[1L]]])

      for (i in seq_along(self$columns)) {
        column_name <- self$columns[i]
        noise <- self$fun(n, self$mean, self$sd)

        s[["result"]][[column_name]] <-
          s[["result"]][[column_name]] + noise

      }
      return(NULL)
    }
  )
)
