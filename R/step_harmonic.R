#' R6 Class
#'
#' `StepHarmonic` generates sin and cosine curves at specified frequencies.
#'
#' @param frequency
#' @param cycle_size
#' @param starting_value
#' @inheritParams Step
#'
#' @export
StepHarmonic <- R6Class(
  classname = 'step_harmonic',
  inherit = Step,

  public = list(

    # step specific variables
    frequency = NA_real_,
    cycle_size = NA_real_,
    starting_value = NA_real_,

    initialize = function(...,
                          frequency = NA_real_,
                          cycle_size = NA_real_,
                          starting_value = NA_real_,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_harmonic"
      type         <- 'add'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)

      # step specific values
      self$frequency <- sort(frequency)
      self$starting_value <- starting_value
      self$cycle_size <- cycle_size

      invisible(self)
    },

    bake = function(new_data) {

      n_frequency <- length(self$frequency)

      column_name     <- self$columns

      hals <- harmonic_list(new_data,
                            frequency = self$frequency,
                            start = self$starting_value,
                            cycle_size = self$cycle_size)

      names(hals) <- file.path(self$id,
                               column_name,
                               rep(c("sin", "cos"), n_frequency),
                               rep(pad_num(n_frequency), each = 2L),
                               fsep = '_')

      hals

    }

  )
)
