#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Adjust the dispersion and central value --------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepNormalize <- R6Class(
  classname = "step_normalize",
  inherit = Step,
  public = list(
    center = c(),
    scale = c(),
    na_rm = TRUE,
    # step specific variables
    initialize = function(terms,
                          role = "predictor",
                          na_rm = TRUE,
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_normalize"
      env_list$type <- "modify"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$na_rm <- na_rm

      invisible(self)
    },
    prep = function(new_data, info) {
      super$prep(new_data, info)

      self$center <- collapse::fmean(unclass(new_data)[self$columns],
        na.rm = self$na_rm,
        drop = TRUE
      )
      self$scale <- collapse::fsd(unclass(new_data)[self$columns],
        na.rm = self$na_rm,
        drop = TRUE
      )
    },
    # subtract the central value from a column
    bake = function(new_data) {

      for (i in seq_along(self$columns)) {
        new_data[[i]] <- (new_data[[i]] - self$center[i]) *
          (1.0 / self$scale[i])
      }

      # fscale(new_data, self$center, self$scale)
      # scale_list_param_std(new_data,
      #                        center = self$center,
      #                        scale = 1.0/self$scale)

      self$result <- new_data
      self$result
    }
  )
)
