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
    # power = NA_real_,
    n_groups = NA_integer_,
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
                          # power = 3,
                          n_groups = 100,
                          time_step = 1.0,
                          formula = NULL,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_fft_transfer_experimental"
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
      # self$power <- power
      self$n_groups <- n_groups
      self$time_step <- time_step


      invisible(self)
    },
    bake = function(new_data, term_info, steps) {

      n  <- length(new_data[[1]])
      df <- 1.0 / n
      frequency <- seq.int(from = 0, by = df, length.out = n) * 86400 / self$time_step
      frequency <- list(frequency = group_frequency(frequency, self$n_groups))
      n_freq <- length(frequency[[1]])

      vars_list <- names(new_data)

      if (!is.null(self$formula)) {
        vars_list <- get_formula_vars(formula = self$formula,
                                      data = unclass(new_data))

      } else {
        vars_list <- list(predictors = self$columns[-1],
                          outcomes = self$columns[1])
      }

      for(i in seq_along(vars_list$outcomes)) {
        tmp_data <- unclass(new_data)[c(vars_list$outcomes[i],
                                        vars_list$predictors)]

        print(c(vars_list$outcomes[i],
                vars_list$predictors))
        print(str(tmp_data))
        res <- collapse::mctl(
          transfer_pgram_smooth(
            collapse::qM(tmp_data),
            self$spans,
            self$detrend,
            self$demean,
            self$taper,
            # self$power,
            self$n_groups
          )
        )

        self$new_columns <- name_columns(self$prefix, NULL, n = length(res))
        names(res) <- self$new_columns
        res <- append(res, frequency)

        res <- append(res, list(variable = rep(vars_list$outcomes[i], n_freq)))
        res <- append(res, list(id = rep(self$id, n_freq)))

        self$fft_result[[i]] <- res

      }

      return(NULL)
    }
  )
)
