#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Calculate the Periodogram ----------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepPgram` calculates the periodogram (estimate of spectral density)
#'
#' @inheritParams Step
#'
#' @export
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
    initialize = function(terms,
                          spans = 3,
                          detrend = TRUE,
                          demean = TRUE,
                          lst = TRUE,
                          taper = 0.1,
                          pad_fft = TRUE,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_fft_pgram"
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
      self$lst <- lst
      self$pad_fft <- pad_fft

      invisible(self)
    },
    bake = function(new_data) {

      if (self$lst) {
        pspec <- collapse::mctl(spec_pgram(
          collapse::qM(new_data),
          self$spans,
          self$detrend,
          self$demean,
          self$taper,
          self$pad_fft
        ))
      } else {
        pspec <- spec_pgram_list(
          new_data,
          self$spans,
          self$detrend,
          self$demean,
          self$taper,
          self$pad_fft
        )
      }
      self$new_columns <- name_columns(self$prefix, NULL, n = length(pspec))
      names(pspec) <- self$new_columns

      return(pspec)
    }
  )
)
