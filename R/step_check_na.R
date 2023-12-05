#' R6 Class
#'
#' `StepCheckNA` adjust the dispersion by standard deviation.
#' @inheritParams Step
#'
#' @export
StepCheckNA <- R6Class(
  classname = 'step_check_na',
  inherit = Step,

  public = list(
    # step specific variables
    initialize = function(...,
                          role = "check",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_check_na"
      type         <- 'check'
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

      chck <- anyNA(new_data)
      names(chck) <- file.path(self$id, self$columns, fsep = "_")
      return(chck)

    }
  )
)

