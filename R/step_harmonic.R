#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Generate Sine and Cosine Terms for Harmonic Analysis -------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
                          starting_value = 0.0,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_harmonic"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
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

        nn <- paste(
          rep(name_columns(self$prefix, column_name, n_frequency), each = 2L),
          rep(c("cos", "sin"), n_frequency),
          sep = "_"
        )

        names(hals[[i]]) <- nn
        self$new_columns <- c(self$new_columns, nn)

      }

      unlist(hals, recursive = FALSE)
    },
    response = function(co) {

      f  <- self$frequency

      nr <- length(f)
      nc <- ncol(co)

      x  <- rep(rep(f, each = nc), 2L)

      cos_coefficient <- co[seq.int(1L, nr * 2L, 2L), , drop = FALSE]
      sin_coefficient <- co[seq.int(2L, nr * 2L, 2L), , drop = FALSE]

      amp_phase <- c(
        as.vector(t(sqrt(cos_coefficient^2 + sin_coefficient^2))), # amplitude
        as.vector(t(atan2(cos_coefficient, sin_coefficient)))      # phase
      )

      variable <- c(
        rep("amplitude", nr * nc),
        rep("phase", nr * nc)
      )

      list(x = x, variable = variable, value = amp_phase, step_id = self$id)
    }

  )
)
