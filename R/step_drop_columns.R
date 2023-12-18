#' R6 Class
#'
#' `StepDropColumns` removes columns from output.
#' @inheritParams Step
#'
#' @export
StepDropColumns <- R6Class(
  classname = 'step_drop_columns',
  inherit = Step,

  public = list(
    initialize = function(...,
                          role = "modify",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_drop_columns"
      type         <- 'modify'
            enq = NULL

      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)


      invisible(self)
    },

    bake = function(new_data) {

      new_data[self$columns] <- list(NULL)

      new_data
    }
  )
)

# a <- sample(1:10, 100000, replace = TRUE)
# bench::mark(funique(a), unique(a))
