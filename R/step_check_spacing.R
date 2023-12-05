#' R6 Class
#'
#' `StepCheckSpacing` adjust the dispersion by standard deviation.
#' @inheritParams Step
#'
#' @export
StepCheckSpacing <- R6Class(
  classname = 'step_check_spacing',
  inherit = Step,

  public = list(
    # step specific variables
    initialize = function(...,
                          role = "check",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_check_spacing"
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

      chck <- collapse::fndistinct(collapse::fdiff(new_data)) == 1L
      names(chck) <- file.path(self$id, self$columns, fsep = "_")

      return(chck)

    }
  )
)

