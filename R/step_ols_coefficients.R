#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Predict Regression Terms -----------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepOlsCoefficients` Uses the Eigen C++ library fast versions to generate
#' predictions from different steps.
#'
#' @inheritParams Step
#'
#' @export
StepOlsCoefficients <- R6Class(
  classname = "step_ols_coefficients",
  inherit = Step,
  public = list(

    # step specific variables
    outcomes = NULL,
    predictors = NULL,
    coefficients = NULL,
    response = NULL,


    initialize = function(terms,
                          role = "augment",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_ols_coefficients"
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

      self$coefficients <- determine_coefficients(x, y)


      return(NULL)
    }
  )
)
