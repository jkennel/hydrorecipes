#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Subset dataset rows using sample ---------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepSubsetNAOmit <- R6Class(
  classname = "step_subset_na_omit",
  inherit = Step,
  public = list(
    # step specific variables
    row_numbers = NULL,
    initialize = function(terms,
                          role = "modify",
                          ...) {

      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_subset_na_omit"
      env_list$type <- "model"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      invisible(self)

    },
    bake = function(r) {

      if (is.null(self$columns)) {
        self$result <- collapse::na_omit(r$get_result())
        r$template_step <- length(r$time_bake) + 1L
        return(NULL)
      }

      res <- r$get_result()
      self$row_numbers <- complete.cases(res[self$columns])
      self$result <- collapse::ss(res, i = self$row_numbers, check = FALSE)
      r$template_step <- length(r$time_bake) + 1L

      return(NULL)
    }
  )
)

