#' R6 Class
#'
#' `StepAddVars` adds variable vectors.
#' @param vars name of vars
#'
#' @inheritParams Step
#'
#' @export
StepAddVars <- R6Class(
  classname = 'step_add_vars',
  inherit = Step,

  public = list(

    # step specific variables

    initialize = function(terms,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- 'step_add_vars'
      env_list$type <- 'add_from_template'
      super$initialize(terms = terms,
                       env_list[names(env_list) != "terms"])

      invisible(self)
    },
    bake = function(new_data) {
      return(unclass(new_data))
    }
  )
)
