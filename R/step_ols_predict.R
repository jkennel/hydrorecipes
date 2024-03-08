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
    response = NULL,
    # residuals = NULL,
    # s = NULL,
    # df_residual = NULL,
    # rank = NULL,
    # std_error = NULL,

    initialize = function(terms,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_ols_predict"
      env_list$type <- "supervised_add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      invisible(self)
    },
    bake = function(new_data, term_info) {

      x <- get_regression_data(new_data, term_info, id_type = "predictor")
      y <- get_regression_data(new_data, term_info, id_type = "outcome")

      self$coefficients <- determine_coefficients(x, y)

      lst <- predict_groups(x, self$coefficients)
      self$new_columns <- name_columns(self$prefix, NULL, length(self$lst))
      names(lst) <- self$new_columns


      # nms <- names(new_data)
      #
      # # term info data
      # ti <- collapse::qDF(term_info)
      # ti <- ti[ti$source != "removed", ]
      # ti <- ti[ti$variable %in% nms, ]
      #
      # # create regression matrices
      # outcomes <- ti[ti$roles == "outcome", ]
      # predictors <- ti[ti$roles == "predictor", ]
      # predictors$inds <- seq_len(nrow(predictors))
      #
      # # save predictor and outcome info
      # self$predictors <- predictors
      # self$outcomes <- outcomes
      #
      # # id data
      # outcome_ids <- which(nms %in% outcomes$variable)
      # predictor_ids <- which(nms %in% predictors$variable)
      #
      # # outcome and predictor data
      # to_rem <- collapse::missing_cases(new_data)
      #
      # m_predictors <- collapse::qM(unclass(new_data)[predictor_ids])
      # m_outcomes <- collapse::qM(unclass(new_data)[outcome_ids])
      #
      # # solve
      # fit <- llt_solve(
      #   m_predictors[!to_rem, , drop = FALSE],
      #   m_outcomes[!to_rem, , drop = FALSE]
      # )
      #
      # self$coefficients <- fit
      #
      # # subsets are the regressor groups
      # subsets <- split(
      #   predictors$inds,
      #   data.table::rleid(predictors$step_index)
      # )
      #
      # lst <- list()
      # for (i in seq_along(subsets)) {
      #   lst[[i]] <- collapse::mctl(
      #     m_predictors[, subsets[[i]], drop = FALSE] %*%
      #       fit[subsets[[i]], , drop = FALSE]
      #   )
      # }
      #
      # names(lst) <- name_columns(self$id, NULL, length(self$lst))
      #
      # lst
    }
  )
)
