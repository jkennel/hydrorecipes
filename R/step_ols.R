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
    coefficients = NULL,
    # decomposition = NULL,
    response_data = NULL,

    do_response = NULL,
    # s = NULL,
    # df_residual = NULL,
    # rank = NULL,
    # std_error = NULL,

    initialize = function(
      formula = NULL,
      role = "predictor",
      do_response = TRUE,
      ...
    ) {
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

      invisible(self)
    },
    bake = function(r) {
      step_names <- r$get_step_field("step_name")
      new_columns <- r$get_step_field("new_columns")

      new_data <- return_type(
        x = r$get_result(),
        type = "m",
        formula = self$formula,
        combined = FALSE
      )

      self$predictors <- new_data[[1L]]
      self$outcomes <- new_data[[2L]]

      co_names <- colnames(self$predictors)
      nms_outcome <- colnames(self$outcomes)
      column_list <- r$get_term_index(co_names)

      to_rem <- !(complete.cases(self$predictors, self$outcomes))

      # steps which have a response
      wh <- which(lengths(column_list) != 0)

      # NEED naming coefficients, fitted.values, decomposition, residuals
      # ols:
      #  - coefficients
      #  - coefficients_list
      #  - fitted.values
      #  - decomposition
      #  - residuals
      #  - s
      #  - df.residual
      #  - rank
      #  - Std. Error

      self$coefficients <- determine_coefficients(
        self$predictors,
        self$outcomes,
        to_rem,
        column_list[wh],
        FALSE
      )$coefficients

      if (self$do_response) {
        self$fit <- predict_decomposition(
          self$predictors,
          self$outcomes,
          self$coefficients,
          subs = column_list[wh]
        )

        resp <- list()
        nms_decomp <- list()

        for (i in seq_along(wh)) {
          co <- self$fit$coefficients_list[[i]]
          colnames(co) <- nms_outcome

          # response
          resp[[i]] <- r$steps[[wh[i]]]$response(co)

          if (resp[[i]]$step_columns[1] == "") {
            nms_decomp[[i]] <- paste(nms_outcome, resp[[i]]$term[1], sep = "_")
          } else {
            nms_decomp[[i]] <- paste(
              nms_outcome,
              paste(resp[[i]]$term[1], resp[[i]]$step_columns[1], sep = "_"),
              sep = "_"
            )
          }
        }

        nms_decomp <- append(
          nms_decomp,
          list(
            paste(nms_outcome, "fitted", sep = "_"),
            paste(nms_outcome, "residuals", sep = "_")
          )
        )

        res <- collapse::rowbind(resp)

        res <- append(res, list(id = rep(self$id, length(res[[1L]]))))
        self$response_data <- res

        names(self$fit$decomposition) <- unlist(nms_decomp)
      }

      return(NULL)
    }
  )
)
