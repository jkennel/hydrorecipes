#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Check a regressor for NA values ----------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepCheckNA`
#' @inheritParams Step
#'
#' @export
StepCheckNA <- R6Class(
  classname = "step_check_na",
  inherit = Step,
  public = list(
    # step specific variables
    initialize = function(terms,
                          role = "check",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_check_na"
      env_list$type <- "check"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"]
      )

      invisible(self)
    },
    bake = function(new_data) {
      chck <- anyNA(new_data)
      names(chck) <- file.path(self$id, self$columns, fsep = "_")
      return(chck)
    }
  )
)
