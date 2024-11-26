# NEED TO FIX SELECTORS TO MAKE EASIER

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Estimate the Coherence Step --------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepCoherence <- R6Class(
  classname = "step_fft_coherence",
  inherit = Step,
  public = list(
    coherence = NULL,

    # step specific variables
    initialize = function(terms,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_fft_coherence"
      env_list$type <- "augment"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      invisible(self)
    },
    bake = function(new_data) {

      self$coherence <- collapse::mctl(ordinary_coherence_phase(collapse::qM(new_data)))
      n <- length(new_data)
      comb <- expand.grid(x = 1:n, y = 1:n)
      comb <- comb[comb$y > comb$x]
      self$new_columns <- paste(self$prefix, comb$x, comb$y, sep = "_")
      setnames(self$coherence, self$new_columns)

      return(NULL)
    }
  )
)
