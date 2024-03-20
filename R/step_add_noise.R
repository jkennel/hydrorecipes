#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Remove the Central Value (mean) from a Regressor Step ------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepAddNoise` add noise to data.
#'
#' @inheritParams Step
#' @inheritParams recipes::step_center
#'
#' @param fun the function to use. Defaults to `rnorm`.
#'
#' @family common
#'
#' @export
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
    bake = function(new_data) {

      for (i in seq_along(self$columns)) {
        noise <- self$fun(self$mean, self$sd)
        new_data[[i]] <- new_data[[i]] + noise
      }

      new_data
    }
  )
)
