#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Remove Regressors Step -------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepDropColumns` removes columns from output.
#' @inheritParams Step
#'
#' @family common
#'
#' @export
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
                       env_list[names(env_list) != "terms"])


      invisible(self)
    },

    bake = function(new_data) {

      columns <- self$columns
      for (i in seq_along(columns)) {
        new_data[columns[i]] <- list(NULL)
      }

      new_data
    }

  )
)
