#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Predict Regression Terms -----------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepOlsPredict` Uses the Eigen C++ library fast versions to generate
#' predictions from different steps.
#'
#' @inheritParams Step
#'
#' @export
StepOlsPredict <- R6Class(
  classname = "step_ols_predict",
  inherit = Step,
  public = list(

    # step specific variables
    # ols_results = list(),
    outcomes = NULL,
    predictors = NULL,
    coefficients = NULL,
    formula = NULL,
    # residuals = NULL,
    # s = NULL,
    # df_residual = NULL,
    # rank = NULL,
    # std_error = NULL,

    initialize = function(formula,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      # terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_ols_predict"
      env_list$type <- "supervise_add"
      super$initialize(
        terms = NULL,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$formula <- formula

      invisible(self)
    },
    bake = function(new_data, term_info) {


      x <- get_regression_data(new_data, term_info, id_type = "predictor")
      y <- get_regression_data(new_data, term_info, id_type = "outcome")


      self$coefficients <- determine_coefficients(x, y)

      # does this work for multiple outcomes?
      lst <- predict_groups(x, self$coefficients)
      lst <- unlist(lst, recursive = FALSE)
      self$new_columns <- paste(names(lst),
                                rep(self$id, length(lst)), sep = "_")
      names(lst) <- self$new_columns
      lst

    }
  )
)
