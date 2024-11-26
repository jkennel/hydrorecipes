#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Dimension Reduction Using Principle Component Analysis -----------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepPca <- R6Class(
  classname = "step_pca",
  inherit = Step,
  public = list(
    pca_results = list(),
    n_comp = NA_integer_,
    na_rm = NA,
    center = NA,
    scale = NA,
    center_values = NA,
    scale_values = NA,

    # step specific variables
    initialize = function(terms,
                          na_rm = TRUE,
                          n_comp = 3,
                          center = TRUE,
                          scale = TRUE,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_pca"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$na_rm <- na_rm
      self$n_comp <- n_comp
      self$center <- center
      self$scale <- scale

      invisible(self)
    },
    prep = function(new_data, info) {
      super$prep(new_data, info)

      new_data <- unclass(new_data)[self$columns]

      if (self$center) {
        self$center_values <- collapse::fmean(new_data, na.rm = self$na_rm)
      } else {
        self$center_values <- rep(0.0, length(new_data))
      }

      if (self$scale) {
        self$scale_values <- collapse::fsd(new_data,
          na.rm = self$na_rm
        )
      } else {
        self$scale_values <- rep(1.0, length(new_data))
      }

      self$pca_results <- pca_list_rotation_eigen(new_data,
        center = self$center_values,
        scale = self$scale_values,
        n_comp = self$n_comp
      )

      invisible(self)

    },
    # subtract the central value from a column
    bake = function(new_data) {

      for (i in seq_along(self$columns)) {
        if (self$center & self$scale) {
          new_data[[i]] <- (new_data[[i]] - self$center_values[i]) *
            (1.0 / self$scale_values[i])
        } else if (self$center) {
          new_data[[i]] <- (new_data[[i]] - self$center_values[i])
        } else if (self$scale) {
          new_data[[i]] <- (new_data[[i]]) * (1.0 / self$scale_values[i])
        }
      }

      new_data <- collapse::qM(new_data)
      new_data <- collapse::mctl(new_data %*% self$pca_results)
      self$new_columns <- name_columns(self$prefix, NULL, n = self$n_comp)
      names(new_data) <- self$new_columns

      self$result <- new_data
      self$result

    }
  )
)
