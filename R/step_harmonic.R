#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Generate Sine and Cosine Terms for Harmonic Analysis -------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepHarmonic` generates sin and cosine curves at specified frequencies.
#'
#' @inheritParams Step
#' @inheritParams recipes::step_harmonic
#'
#' @export
StepHarmonic <- R6Class(
  classname = "step_harmonic",
  inherit = Step,
  public = list(

    # step specific variables
    frequency = NA_real_,
    cycle_size = NA_real_,
    starting_value = NA_real_,
    initialize = function(terms,
                          frequency = NA_real_,
                          cycle_size = NA_real_,
                          starting_value = NA_real_,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_harmonic"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"]
      )


      # step specific values
      self$frequency <- sort(frequency)
      self$starting_value <- starting_value
      self$cycle_size <- cycle_size

      invisible(self)
    },
    bake = function(new_data) {
      n_frequency <- length(self$frequency)

      column_name <- self$columns

      hals <- list()
      for (i in seq_along(column_name)) {
        hals[[i]] <- harmonic_list(unclass(new_data)[[i]],
          frequency = self$frequency,
          start = self$starting_value,
          cycle_size = self$cycle_size
        )
        nms <- name_columns(self$id, column_name, n_frequency)
        names(hals[[i]]) <- paste(rep(nms, each = 2L),
          rep(c("sin", "cos"), n_frequency),
          sep = "_"
        )
      }

      unlist(hals, recursive = FALSE)
    }
  )
)
