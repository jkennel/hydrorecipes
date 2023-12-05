#' R6 Class
#'
#' `StepVarying` adjust the central value to zero.
#' @inheritParams Step
#'
#' @export
StepVarying <- R6Class(
  classname = 'step_varying',
  inherit = Step,

  public = list(
    to_remove = NULL,

    # step specific variables
    initialize = function(...,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_varying"
      type         <- 'modify'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)

      invisible(self)
    },

    # subtract the central value from a column
    bake = function(new_data) {

      self$to_remove <- !varying(new_data)
      new_data[self$to_remove] <- list(NULL)

      new_data
    }
  )
)

