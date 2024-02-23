#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Add an Intercept Term --------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepIntercept` adds variable vectors.
#'
#' @inheritParams Step
#'
#' @export
StepIntercept <- R6Class(
  classname = "step_intercept",
  inherit = Step,
  public = list(

    # step specific variables
    initialize = function(terms,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      env_list <- get_function_arguments()
      env_list$step_name <- "step_intercept"
      env_list$type <- "add"
      super$initialize(
        terms = NULL,
        env_list[names(env_list) != "terms"]
      )


      invisible(self)
    },
    bake = function(new_data) {
      return(setNames(list(rep(1.0, length(new_data[[1]]))), self$id))
    }
  )
)
