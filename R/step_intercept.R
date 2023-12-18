#' R6 Class
#'
#' `StepIntercept` adds variable vectors.
#' @param vars name of vars
#'
#' @inheritParams Step
#'
#' @export
StepIntercept <- R6Class(
  classname = 'step_intercept',
  inherit = Step,

  public = list(

    # step specific variables

    initialize = function(...,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_intercept"
      type         <- 'add'
      enq          <- NULL
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      print(str(inputs))
      do.call(super$initialize, inputs)

      invisible(self)
    },
    bake = function(new_data) {
      return(list(intercept = rep(1.0, length(new_data))))
    }

  )
)
