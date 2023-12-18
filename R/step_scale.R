#' R6 Class
#'
#' `StepScale` adjust the dispersion by the standard deviation.
#'
#' @inheritParams Step
#' @inheritParams recipes::step_scale
#' @param fun the function to use for calculating the distpersion. The default is
#' `collapse::fsd`
#'
#' @export
StepScale <- R6Class(
  classname = 'step_scale',
  inherit = Step,

  public = list(
    column_values = c(),
    na_rm = NA,
    fun = NULL,
    n_sd = NA_integer_,

    # step specific variables
    initialize = function(...,
                          role = "predictor",
                          skip = FALSE,
                          na_rm = TRUE,
                          fun = collapse::fsd,
                          n_sd = 1L,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_scale"
      type         <- 'modify'
      enq          <- NULL
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)

      self$na_rm <- na_rm
      self$fun <- fun
      self$n_sd <- n_sd

      invisible(self)
    },

    prep = function(new_data) {
      self$column_values <- self$fun(unclass(new_data)[self$columns],
                                     na.rm = self$na_rm, drop = TRUE) * self$n_sd

      self$column_values <- 1.0 / self$column_values
    },

    # subtract the central value from a column
    bake = function(new_data) {

      for(i in seq_along(self$columns)) {
        new_data[[i]] %*=% self$column_values[i]
      }

      new_data
    }
  )
)

