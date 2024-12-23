#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate Welch's Periodogram  -----------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepWelch <- R6Class(
  classname = "step_fft_welch",
  inherit = Step,
  public = list(

    # step specific variables
    length_subset = NA,
    overlap = NA_real_,
    window = NULL,
    time_step = NA_real_,
    fft_result = NA,
    initialize = function(terms,
                          length_subset,
                          overlap = 0.8,
                          window,
                          time_step = 1.0,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_fft_welch"
      env_list$type <- "augment"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$length_subset <- length_subset
      self$overlap <- overlap
      self$window <- window
      self$time_step <- time_step

      invisible(self)
    },
    bake = function(s) {
      self$fft_result <- collapse::mctl(spec_welch(
        collapse::qM(s[["result"]][self$columns]),
        self$length_subset,
        self$overlap,
        self$window
      ))

      self$new_columns <- name_columns(self$prefix, "", n = length(self$fft_result))

      names(self$fft_result) <- self$new_columns

      n  <- length(self$fft_result[[1L]])
      df <- 1.0 / n
      frequency <- list(frequency = seq.int(from = 0.0, by = df,
                                            length.out = n) * 86400.0 / self$time_step)
      self$fft_result <- append(self$fft_result, frequency)


      return(NULL)
    }
  )
)
