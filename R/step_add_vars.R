#' R6 Class
#'
#' `StepAddVars` adds variable vectors.
#' @param vars name of vars
#' @inheritParams Step
#'
#' @export
StepAddVars <- R6Class(
  classname = 'step_add_vars',
  inherit = Step,

  public = list(
    vars = NULL,
    # step specific variables
    initialize = function(...,
                          vars,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_add_vars"
      type         <- 'add'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)

      self$vars <- vars

      invisible(self)
    },
    bake = function(new_data) {
      return(unclass(new_data)[unique(self$vars)])
    }
  )
)
