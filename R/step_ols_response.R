#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Predict Regression Terms -----------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepOlsResponse` Uses the Eigen C++ library fast versions to generate
#' predictions from different steps.
#'
#' @inheritParams Step
#'
#' @export
StepOlsResponse <- R6Class(
  classname = "step_ols_response",
  inherit = Step,
  public = list(

    # step specific variables
    # ols_results = list(),
    outcomes = NULL,
    predictors = NULL,
    coefficients = NULL,
    response_data = NULL,
    # residuals = NULL,
    # s = NULL,
    # df_residual = NULL,
    # rank = NULL,
    # std_error = NULL,

    initialize = function(terms,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_ols_response"
      env_list$type <- "supervise_augment"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      invisible(self)
    },
    bake = function(new_data, term_info, steps) {
      x <- get_regression_data(new_data, term_info, id_type = "predictor")
      y <- get_regression_data(new_data, term_info, id_type = "outcome")
      print(str(x))
      print(str(y))

      self$coefficients <- determine_coefficients(x, y)

      self$response_data <- response_groups(steps, x, self$coefficients)

      return(NULL)
    }
  )
)
