#' R6 Class
#'
#' `StepDummy` dummy encoding for factor or integer input.
#'
#' @inheritParams Step
#'
#' @export
StepDummy <- R6Class(
  classname = 'step_dummy',
  inherit = Step,

  public = list(

    # step specific variables
    levels = NULL,

    #' @description
    #' @inheritParams StepAddVars
    #' @return A new `Step`.
    initialize = function(...,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_dummy"
      type         <- 'add'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)

      invisible(self)
    },
    prep = function(new_data) {
      self$levels <- levels(new_data)

    },
    bake = function(new_data) {

      column_name <- self$columns

      # check for new level(s)
      if(sum(levels(new_data) %!in% self$levels) > 0L) {
        warning(file.path("New levels found during bake step. (", self$id, ")", fsep = ""))
      }

      dum <- to_dummy(new_data)

      names(dum) <- file.path(self$id,
                              column_name,
                              pad_num(length(dum)),
                              fsep = "_")

      dum
    }

  )
)

