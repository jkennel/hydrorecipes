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
      env_list$type <- "model"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      invisible(self)
    },
    bake = function(r) {

      tf <- r$get_transfer_data()
      nms <- names(tf)

      self$coherence <- lapply(tf, function(x) {

        nms <- names(x)
        x <- collapse::qM(x[-which(nms == "frequency")])
        res <- collapse::mctl(ordinary_coherence_phase(x))
        n <- get_column_number(length(res))

        comb <- expand.grid(x = 1:n, y = 1:n)
        comb <- comb[comb$y > comb$x, ]

        self$new_columns <- paste(self$prefix,
                                  c(rep("coherence", n),
                                  rep("phase", n)),
                                  comb$x, comb$y, sep = "_")
        # print(str(res))
        # print(self$new_columns)
        res <- setNames(res, self$new_columns)
        res

      })


      return(NULL)
    }
  )
)
