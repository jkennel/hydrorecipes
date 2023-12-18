#' R6 Class
#'
#' `StepPredictOLS` Uses the Eigen C++ library fast versions to generate
#' predictions from different steps.
#'
#' @inheritParams Step
#'
#' @export
StepPredictOLS <- R6Class(
  classname = 'step_predict_ols',
  inherit = Step,

  public = list(

    # step specific variables
    # ols_results = list(),
    outcomes = NULL,
    predictors = NULL,
    coefficients = NULL,
    # residuals = NULL,
    # s = NULL,
    # df_residual = NULL,
    # rank = NULL,
    # std_error = NULL,

    initialize = function(...,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_predict_ols"
      type         <- 'supervised_add'
      enq <- NULL
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)
      invisible(self)
    },
    # subtract the central value from a column
    bake = function(new_data, term_info) {


      # remove na values
      nms <- names(new_data)


      # term info data
      ti  <- qDF(term_info)
      ti  <- ti[ti$variable %in% nms,]


      # create regression matrices
      outcomes   <- ti[ti$roles == "outcome", ]
      predictors <- ti[ti$roles == "predictor", ]
      predictors$inds <- 1:nrow(predictors)
      subsets <- split(predictors$inds, data.table::rleid(predictors$step_index))

      # save predictor and outcome info
      self$predictors <- predictors
      self$outcomes <- outcomes

      # id data
      outcome_ids <- which(nms %in% outcomes$variable)
      predictor_ids <- which(nms %in% predictors$variable)


      # outcome and predictor data
      to_rem <- missing_cases(new_data)
      # no_na <- new_data[!to_rem,,drop = FALSE]
      m_predictors <- collapse::qM(new_data[predictor_ids])
      m_outcomes <- collapse::qM(new_data[outcome_ids])

      # solve
      fit <- llt_solve(m_predictors[!to_rem, , drop = FALSE],
                       m_outcomes[!to_rem, , drop = FALSE])


      self$coefficients <- fit


      lst <- list()
      for (i in seq_along(subsets)) {

        lst[[i]] <- m_predictors[, subsets[[i]], drop = FALSE] %*%
          fit[subsets[[i]], , drop = FALSE]

      }

      names(lst) <- file.path(self$id,
                              pad_num(length(lst)),
                              fsep = "_")

      lst

    }
  )
)

