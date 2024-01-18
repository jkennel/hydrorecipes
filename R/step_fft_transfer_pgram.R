#' R6 Class
#'
#' `StepTransferPgram` adds variable vectors.
#' @param vars name of vars
#'
#' @inheritParams Step
#'
#' @export
StepTransferPgram <- R6Class(
  classname = 'step_fft_transfer_pgram',
  inherit = Step,

  public = list(

    # step specific variables
    spans = NA_integer_,
    detrend = NA,
    demean = NA,
    taper = NA_real_,
    initialize = function(terms,
                          spans = 3,
                          detrend = TRUE,
                          demean = TRUE,
                          taper = 0.1,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- 'step_fft_transfer_pgram'
      env_list$type <- 'add'
      super$initialize(terms = terms,
                       env_list[names(env_list) != "terms"])

      self$spans <- spans
      self$detrend <- detrend
      self$demean <- demean
      self$taper <- taper

      invisible(self)
    },
    bake = function(new_data) {
      tf <- collapse::mctl(
        transfer_pgram(collapse::qM(new_data),
                       self$spans,
                       self$detrend,
                       self$demean,
                       self$taper))
    }

  )
)
