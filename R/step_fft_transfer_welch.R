#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate the Transfer Function from Periodograms  ---------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepTransferWelch <- R6Class(
  classname = "step_fft_transfer_welch",
  inherit = Step,
  public = list(

    # step specific variables
    length_subset = NA,
    overlap = NA_real_,
    window = NA_real_,
    time_step = NA_real_,

    formula = NULL,
    outcomes = NULL,
    predictors = NULL,

    fft_result = list(),

    initialize = function(terms,
                          length_subset,
                          overlap = 0.8,
                          window,
                          time_step = 1.0,
                          formula = NULL,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_fft_transfer_welch"
      env_list$type <- "supervise_augment"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$formula <- formula
      self$length_subset <- length_subset
      self$overlap <- overlap
      self$window <- window
      self$time_step <- time_step

      invisible(self)
    },
    bake = function(new_data, term_info, steps) {

      vars_list <- select_fft_vars_list(new_data, self$formula, self$columns)

      for(i in seq_along(vars_list$outcomes)) {

        tmp_data <- unclass(new_data)[c(vars_list$outcomes[i],
                                        vars_list$predictors)]

        res <- collapse::mctl(
          transfer_welch(
            collapse::qM(tmp_data),
            self$length_subset,
            self$overlap,
            self$window
          )
        )

        n  <- length(res[[1]])
        df <- 1 / n
        frequency <- list(frequency = seq.int(from = 0, by = df,
                                              length.out = n) * 86400 / self$time_step)

        self$new_columns <- name_columns(paste(names(tmp_data), collapse = "_"), NULL, n = length(res))

        names(res) <- self$new_columns
        res <- append(res, frequency)

        # res <- append(res, list(variable = rep(vars_list$outcomes[i], n)))
        res <- append(res, list(id = rep(self$id, n)))

        self$fft_result[[i]] <- res
      }

      return(NULL)
    }
  )
)
