#' R6 Class
#'
#' `StepCenter` adjust the central value to zero.
#' @inheritParams Step
#'
#' @export
StepCenter <- R6Class(
  classname = 'step_center',
  inherit = Step,

  public = list(
    column_values = c(),
    na_rm = NA,
    fun = NULL,

    # step specific variables
    initialize = function(...,
                          role = "predictor",
                          skip = FALSE,
                          na_rm = TRUE,
                          fun = collapse::fmean,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_center"
      type         <- 'modify'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)

      self$na_rm <- na_rm
      self$fun <- fun

      invisible(self)
    },

    prep = function(new_data) {
      self$column_values <- self$fun(unclass(new_data)[self$columns],
                                     na.rm = self$na_rm, drop = TRUE)
    },

    # subtract the central value from a column
    bake = function(new_data) {

      for(i in seq_along(self$columns)) {
        new_data[[i]] %-=% self$column_values[i]
      }

      new_data
    }
  )
)

