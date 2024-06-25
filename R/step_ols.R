#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Predict Regression Terms -----------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepOls <- R6Class(
  classname = "step_ols",
  inherit = Step,
  public = list(

    # step specific variables
    outcomes = NULL,
    predictors = NULL,
    coefficients = NULL,
    formula = NULL,
    decomposition = NULL,
    response_data = NULL,

    do_response = NULL,
    do_predict = NULL,
    # s = NULL,
    # df_residual = NULL,
    # rank = NULL,
    # std_error = NULL,

    initialize = function(formula,
                          role = "predictor",
                          do_response = TRUE,
                          do_predict = TRUE,
                          ...) {
      # get function parameters to pass to parent
      # terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_ols"
      env_list$type <- "supervise_augment"
      super$initialize(
        terms = NULL,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$formula <- formula
      self$do_response <- do_response
      self$do_predict <- do_predict

      invisible(self)
    },
    bake = function(new_data, term_info, steps) {


      self$predictors <- get_regression_data(new_data, term_info, id_type = "predictor")
      self$outcomes <- get_regression_data(new_data, term_info, id_type = "outcome")

      self$coefficients <- determine_coefficients(self$predictors, self$outcomes)

      if (self$do_predict) {

        # predict for each group
        self$decomposition <- predict_groups(self$predictors, self$coefficients, steps)
        self$decomposition <- unlist(self$decomposition, recursive = FALSE)
        self$decomposition <- append(self$decomposition,
                                     list(id = rep(self$id, length(self$decomposition[[1]]))))
      }

      if (self$do_response) {
        # response for each group
        # column names in term info
        co_names <- self$predictors$term_info$variable

        resp <- list()

        for (i in seq_along(steps)) {

          wh  <- collapse::whichv(self$predictors$term_info$step_index, i)
          co_name <- co_names[wh]

          if (length(co_name) > 0) {
            co <- self$coefficients[wh, , drop = FALSE]
            resp[[i]] <- steps[[i]]$response(co)
            if (!"outcome" %in% names(resp[[i]])) {
              resp[[i]]$outcome <- rep(colnames(co), times = nrow(co))
            }
            if (!"term" %in% names(resp[[i]])) {
              resp[[i]]$term <- rep(co_name, times = ncol(co))
            }
            resp[[i]]$step_columns <- paste(steps[[i]]$columns, collapse = "_")
          }
        }


        # save the response
        self$response_data <- collapse::rowbind(resp)
        self$response_data <- append(self$response_data,
                                     list(id = rep(self$id, length(self$response_data[[1]]))))
      }

      return(NULL)


    }
  )
)
