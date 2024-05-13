#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate the Transfer Function from Periodograms  ---------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepTransferWelch <- R6Class(
  classname = "step_transfer_welch",
  inherit = Step,
  public = list(

    # step specific variables
    length_subset = NA,
    overlap = NA_real_,
    window = NA_real_,
    initialize = function(terms,
                          length_subset,
                          overlap = 0.8,
                          window,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_transfer_welch"
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
      tf <- collapse::mctl(
        transfer_welch(
          collapse::qM(new_data),
          self$length_subset,
          self$overlap,
          self$window
        )
      )
      self$new_columns <- name_columns(self$prefix, NULL, n = length(tf))

      names(tf) <- self$new_columns
    }
  )
)
