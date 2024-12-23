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
    fit = NULL,
    formula = NULL,
    # decomposition = NULL,
    response_data = NULL,

    do_response = NULL,
    # do_predict = NULL,
    # s = NULL,
    # df_residual = NULL,
    # rank = NULL,
    # std_error = NULL,

    initialize = function(formula = NULL,
                          role = "predictor",
                          do_response = TRUE,
                          # do_predict = TRUE,
                          ...) {
      # get function parameters to pass to parent
      # terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_ols"
      env_list$type <- "model"
      super$initialize(
        terms = NULL,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$formula <- formula
      self$do_response <- do_response
      # self$do_predict <- do_predict

      invisible(self)
    },
    bake = function(r) {

      step_names  <- r$get_step_field("step_name")

      new_data <- return_type(x = r$get_result(),
                              type = "m",
                              formula = self$formula,
                              combined = FALSE)

      self$predictors <- new_data[[1L]]
      self$outcomes   <- new_data[[2L]]

      co_names <- colnames(self$predictors)
      nms_outcome <- colnames(self$outcomes)
      column_list <- r$get_term_index(co_names)

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
                                         column_list)



      if (self$do_response) {
        new_columns <- r$get_step_field("new_columns")
        # new_columns <- new_columns[lengths(new_columns) != 0]

        # response for each group
        # column names in term info

        resp <- list()

        # print("xxxxxxxx")
        # print(new_columns)
        for (i in seq_along(new_columns)) {
          nc <- new_columns[[i]]
          if (is.null(nc)) {
            next
          }

          wh  <- co_names %in% nc
          if(sum(wh) == 0) {
            next
          }

          co_name <- co_names[wh]
          # print("---------")
          # print(wh)
          # print(co_name)
          if (length(co_name) > 0) {

            co <- self$fit$coefficients[wh, , drop = FALSE]
            # print(co)
            colnames(co) <- nms_outcome
            resp[[i]] <- r$steps[[i]]$response(co)
            # print(str(resp[[i]]))

            if (!"outcome" %in% names(resp[[i]])) {
              resp[[i]]$outcome <- rep(nms_outcome, times = nrow(co))
            }

            if (!"term" %in% names(resp[[i]])) {
              resp[[i]]$term <- rep(co_name, times = ncol(co))
            }

            resp[[i]]$step_columns <- paste(r$steps[[i]]$columns, collapse = "_")
          }
        }

        res <- collapse::rowbind(resp)

        res <- append(res, list(id = rep(self$id, length(res[[1L]]))))
        self$response_data <- res

      }

      return(NULL)


    }
  )
)
