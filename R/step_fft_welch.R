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
    initialize = function(terms,
                          length_subset,
                          overlap = 0.8,
                          window,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_fft_welch"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$length_subset <- length_subset
      self$overlap <- overlap
      self$window <- window

      invisible(self)
    },
    bake = function(new_data) {
      pspec <- collapse::mctl(spec_welch(
        collapse::qM(new_data),
        self$length_subset,
        self$overlap,
        self$window
      ))

      self$new_columns <- name_columns(self$prefix, "", n = length(pspec))

      names(pspec) <- self$new_columns

      return(pspec)
    }
  )
)
