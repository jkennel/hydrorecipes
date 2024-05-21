#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate the Transfer Function from Periodograms  ---------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepTransferExperimental <- R6Class(
  classname = "step_fft_transfer_experimental",
  inherit = Step,
  public = list(

    # step specific variables
    spans = NA_integer_,
    detrend = NA,
    demean = NA,
    taper = NA_real_,
    power = NA_real_,
    n_groups = NA_integer_,
    time_step = NA_real_,

    fft_result = NA,

    initialize = function(terms,
                          spans = 3,
                          detrend = TRUE,
                          demean = TRUE,
                          taper = 0.1,
                          power = 3,
                          n_groups = 200,
                          time_step = 1.0,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_fft_transfer_experimental"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$spans <- spans
      self$detrend <- detrend
      self$demean <- demean
      self$taper <- taper
      self$power <- power
      self$n_groups <- n_groups
      self$time_step <- time_step


      invisible(self)
    },
    bake = function(new_data) {
      self$fft_result <- collapse::mctl(
        transfer_pgram_smooth(
          collapse::qM(new_data),
          self$spans,
          self$detrend,
          self$demean,
          self$taper,
          self$power,
          self$n_groups
        )
      )
      self$new_columns <- name_columns(self$prefix, NULL, n = length(self$fft_result))
      names(self$fft_result) <- self$new_columns

      n  <- length(new_data[[1]])
      df <- 1.0 / n
      frequency <- seq.int(from = 0, by = df,
                                            length.out = n) * 86400 / self$time_step
      # print(head(frequency))
      # print((frequency))
      # print(hydrorecipes:::group_frequency(frequency, 200))

      frequency <- list(frequency = group_frequency(frequency, self$n_groups))
      self$fft_result <- append(self$fft_result, frequency)

      return(NULL)
    }
  )
)
