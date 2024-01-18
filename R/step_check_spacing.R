#' R6 Class
#'
#' `StepCheckSpacing`
#' @inheritParams Step
#'
#' @export
StepCheckSpacing <- R6Class(
  classname = 'step_check_spacing',
  inherit = Step,

  public = list(
    # step specific variables
    initialize = function(terms,
                          role = "check",
                          ...) {

     # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_check_spacing"
      env_list$type <- "check"
      super$initialize(terms = terms,
                       env_list[names(env_list) != "terms"])

      invisible(self)
    },

    bake = function(new_data) {

      chck <- collapse::fndistinct(collapse::fdiff(new_data)) == 1L
      names(chck) <- file.path(self$id, self$columns, fsep = "_")

      return(chck)

    }
  )
)

