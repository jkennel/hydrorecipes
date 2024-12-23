#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate the Transfer Function from Periodograms  ---------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepTransferPgram <- R6Class(
  classname = "step_fft_transfer_pgram",
  inherit = Step,
  public = list(

    # step specific variables
    spans = NA_integer_,
    detrend = NA,
    demean = NA,
    taper = NA_real_,
    time_step = NA_real_,

    formula = NULL,
    outcomes = NULL,
    predictors = NULL,

    fft_result = list(),

    initialize = function(terms,
                          spans = 3,
                          detrend = TRUE,
                          demean = TRUE,
                          taper = 0.1,
                          time_step = 1.0,
                          formula = NULL,
                          role = "augment",
                          ...) {

      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_fft_transfer_pgram"
      env_list$type <- "supervise_augment"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )


      self$formula <- formula
      self$spans <- spans
      self$detrend <- detrend
      self$demean <- demean
      self$taper <- taper
      self$time_step <- time_step


      invisible(self)
    },
    bake = function(s) {

      vars_list <- select_fft_vars_list(s[["result"]], self$formula, self$columns)

      for(i in seq_along(vars_list$outcomes)) {

        tmp_data <- s[["result"]][c(vars_list$outcomes[i],
                                    vars_list$predictors)]

        res <- collapse::mctl(
          transfer_pgram(
            collapse::qM(tmp_data),
            self$spans,
            self$detrend,
            self$demean,
            self$taper
          )
        )

        n  <- length(res[[1L]])
        df <- 1.0 / n
        frequency <- list(frequency = seq.int(from = 0.0, by = df,
                                              length.out = n) * 86400.0 / self$time_step)

        self$new_columns <- name_columns(paste(names(tmp_data), collapse = "_"), NULL, n = length(res))
        names(res) <- self$new_columns
        res <- append(res, frequency)

        # res <- append(res, list(variable = rep(vars_list$outcomes[i], n)))
        res <- append(res, list(id = rep.int(self$id, n)))

        self$fft_result[[i]] <- res

      }

      return(NULL)
    }
  )
)
