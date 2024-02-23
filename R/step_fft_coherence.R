# NEED TO FIX SELECTORS

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Estimate the Coherence Step --------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepCoherence` estimates the coherence between terms.
#'
#' @inheritParams Step
#'
#' @export
StepCoherence <- R6Class(
  classname = "step_fft_coherence",
  inherit = Step,
  public = list(

    # step specific variables
    initialize = function(terms,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_fft_coherence"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"]
      )

      invisible(self)
    },
    bake = function(new_data) {
      print(new_data)
      collapse::mctl(ordinary_coherence_phase(collapse::qM(new_data)))
    }
  )
)
