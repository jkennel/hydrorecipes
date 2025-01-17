#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Subset dataset rows using sample ---------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepSubsetSample <- R6Class(
  classname = "step_subset_sample",
  inherit = Step,
  public = list(
    # step specific variables
    size = NULL,
    row_numbers = NULL,
    initialize = function(terms,
                          size = 1L,
                          role = "modify",
                          ...) {

      # get function parameters to pass to parent
      env_list <- get_function_arguments()
      env_list$step_name <- "step_subset_sample"
      env_list$type <- "model"

      super$initialize(
        terms = NULL,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$size <- as.integer(size)

      invisible(self)

    },
    bake = function(r) {

      n <- length(r$steps[[r$template_step]]$result[[1L]])
      self$row_numbers <- sample(1:n, self$size, replace = TRUE)

      self$result <- collapse::ss(r$get_result(type = "list"),
                                  i = self$row_numbers, check = FALSE)

      r$template_step <- length(r$time_bake) + 1L

      return(NULL)
    }
  )
)

