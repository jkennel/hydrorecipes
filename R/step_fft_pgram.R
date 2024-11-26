#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate the Periodogram ----------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepPgram <- R6Class(
  classname = "step_fft_pgram",
  inherit = Step,
  public = list(

    # step specific variables
    spans = NA_integer_,
    detrend = NA,
    demean = NA,
    taper = NA_real_,
    lst = NA,
    pad_fft = NA,
    time_step = NA,

    fft_result = list(),

    initialize = function(terms,
                          spans = 3,
                          detrend = TRUE,
                          demean = TRUE,
                          lst = TRUE,
                          taper = 0.1,
                          pad_fft = TRUE,
                          time_step = 1,
                          role = "augment",
                          ...) {

      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_fft_pgram"
      env_list$type <- "augment"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$spans <- spans
      self$detrend <- detrend
      self$demean <- demean
      self$taper <- taper
      self$lst <- lst
      self$pad_fft <- pad_fft
      self$time_step <- time_step

      invisible(self)
    },
    bake = function(new_data) {

      if (self$lst) {
        self$fft_result <- collapse::mctl(spec_pgram(
          collapse::qM(new_data),
          self$spans,
          self$detrend,
          self$demean,
          self$taper,
          self$pad_fft
        ))
      } else {
        self$fft_result <- spec_pgram_list(
          new_data,
          self$spans,
          self$detrend,
          self$demean,
          self$taper,
          self$pad_fft
        )
      }

      self$new_columns <- name_columns(self$prefix, NULL, n = length(self$fft_result))
      names(self$fft_result) <- self$new_columns


      n  <- length(self$fft_result[[1]])
      df <- 1 / n
      frequency <- list(frequency = seq.int(from = 0, by = df,
                           length.out = n) * 86400 / self$time_step)
      # print(str(self$fft_result))
      # print(str(frequency))
      self$fft_result <- append(self$fft_result, frequency)
      # print(is.list(self$fft_result))

      return(NULL)
    }
  )
)
