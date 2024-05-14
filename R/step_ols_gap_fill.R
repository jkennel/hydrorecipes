#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Fill in Gaps using Regression ------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepOlsGapFill <- R6Class(
  classname = "step_ols_gap_fill",
  inherit = Step,
  public = list(

    # step specific variables
    recipe = NULL,
    coefficients = NULL,
    initialize = function(terms,
                          recipe,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_ols_gap_fill"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )


      # step specific values
      self$recipe <- recipe

      invisible(self)
    },
    bake = function(new_data) {

      rec <- self$recipe
      rec <- rec$prep()$bake(data = new_data)
      dat <- rec$result

      x <- get_regression_data(dat, rec$term_info, id_type = "predictor")
      y <- get_regression_data(dat, rec$term_info, id_type = "outcome")

      mode(x$data) <- "double"
      mode(y$data) <- "double"


      self$coefficients <- determine_coefficients(x, y)

      print(str(x$data))
      lst <- collapse::mctl(x$data[, , drop = FALSE] %*% self$coefficients[, , drop = FALSE])

      self$new_columns <- name_columns(self$prefix, colnames(y$data), n = ncol(y$data))
      names(lst) <- self$new_columns

      lst

    }
  )
)



