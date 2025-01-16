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
    fit = NULL,
    predictors = NULL,
    outcomes = NULL,

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
    bake = function(s) {

      r <- self$recipe

      new_data <- return_type(x = r$get_result(),
                              type = "m",
                              formula = r$formula,
                              combined = FALSE)

      self$predictors <- new_data[[1L]]
      self$outcomes   <- new_data[[2L]]

      co_names <- colnames(self$predictors)
      nms_outcome <- colnames(self$outcomes)
      column_list <- r$get_term_index(co_names)
      wh <- which(lengths(column_list) != 0)

      to_rem <- !(complete.cases(self$predictors, self$outcomes))


      # NEED naming coefficients, fitted.values, decomposition, residuals
      # ols:
      #  - coefficients
      #  - fitted.values
      #  - decomposition
      #  - residuals
      #  - s
      #  - df.residual
      #  - rank
      #  - Std. Error
      self$fit <- determine_coefficients(self$predictors,
                                         self$outcomes,
                                         to_rem,
                                         column_list[wh],
                                         FALSE)

      lst <- collapse::mctl(self$predictors %*% self$fit$coefficients)

      self$new_columns <- name_columns(self$prefix, colnames(self$outcomes), n = ncol(self$outcomes))
      names(lst) <- self$new_columns

      self$result <- lst

      return(NULL)

    }
  )
)



