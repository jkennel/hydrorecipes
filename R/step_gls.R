#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Predict Regression Terms -----------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepGls <- R6Class(
  classname = "step_gls",
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

    family = c(
      "gaussian",
      "binomial",
      "poisson",
      "multinomial",
      "cox",
      "mgaussian"
    ),
    weights = NULL,
    offset = NULL,
    alpha = 1,
    nlambda = 100L,
    lambda.min.ratio = NULL,
    lambda = NULL,
    standardize = TRUE,
    intercept = TRUE,
    thresh = 1.0e-7,
    dfmax = NULL,
    pmax = NULL,
    exclude = NULL,
    penalty.factor = NULL,
    lower.limits = -Inf,
    upper.limits = Inf,
    maxit = 1e5,
    type.gaussian = NULL,
    type.logistic = c("Newton", "modified.Newton"),
    standardize.response = FALSE,
    type.multinomial = c("ungrouped", "grouped"),
    relax = FALSE,
    trace.it = 0,
    cox.ties = c("breslow", "efron"),
    control = list(),
    s = NULL,
    # df_residual = NULL,
    # rank = NULL,
    # std_error = NULL,

    initialize = function(
      formula = NULL,
      role = "predictor",
      do_response = TRUE,
      family = c(
        "gaussian",
        "binomial",
        "poisson",
        "multinomial",
        "cox",
        "mgaussian"
      ),
      weights = NULL,
      offset = NULL,
      alpha = 1,
      nlambda = 100L,
      lambda.min.ratio = NULL,
      lambda = NULL,
      standardize = TRUE,
      intercept = FALSE,
      # thresh = 1.0e-7,
      # dfmax = NULL,
      # pmax = NULL,
      exclude = NULL,
      penalty.factor = NULL,
      lower.limits = -Inf,
      upper.limits = Inf,
      # maxit = 1e5,
      type.gaussian = "covariance",
      type.logistic = c("Newton", "modified.Newton"),
      standardize.response = FALSE,
      type.multinomial = c("ungrouped", "grouped"),
      relax = FALSE,
      # trace.it = 0,
      cox.ties = c("breslow", "efron"),
      control = list(),
      s = 0.0,
      ...
    ) {
      # get function parameters to pass to parent
      # terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_gls"
      env_list$type <- "model"
      super$initialize(
        terms = NULL,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$formula <- formula
      self$do_response <- do_response

      self$family = family
      self$offset = offset
      self$alpha = alpha
      self$nlambda = nlambda
      self$lambda.min.ratio = lambda.min.ratio
      self$lambda = lambda
      self$standardize = standardize
      self$intercept = intercept
      # self$thresh = thresh
      # self$dfmax = dfmax
      # self$pmax = pmax
      self$exclude = exclude
      self$penalty.factor = penalty.factor
      self$lower.limits = lower.limits
      self$upper.limits = upper.limits
      # self$maxit = maxit
      self$type.gaussian = type.gaussian
      self$type.logistic = type.logistic
      self$standardize.response = standardize.response
      self$type.multinomial = type.multinomial
      self$relax = relax
      # self$trace.it = trace.it
      self$cox.ties = cox.ties
      self$control = control
      self$s <- s

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
      n <- ncol(self$predictors)
      inter <- "intercept" %in% colnames(self$predictors)
      n <- n - inter

      if (is.null(self$penalty.factor)) {
        self$penalty.factor <- rep(1.0, n)
      }
      if (is.null(self$lambda.min.ratio)) {
        self$lambda.min.ratio <- 1e-4
      }
      if (is.null(self$type.gaussian)) {
        self$type.gaussian <- ifelse(n < 500, "covariance", "naive")
      }
      # need to pass alpha and other arguments !!!!

      self$coefficients <- determine_coefficients_gls(
        self$predictors,
        self$outcomes,
        to_rem,
        column_list[wh],
        FALSE,

        weights = self$weights,
        offset = self$offset,
        alpha = self$alpha,
        nlambda = self$nlambda,
        nlambda.min.ratio = self$nlambda.min.ratio,
        lambda = self$lambda,
        standardize = self$standardize,
        intercept = self$intercept,
        # thresh = self$thresh,
        # pmax = self$pmax,
        exclude = self$exclude,
        penalty.factor = self$penalty.factor,
        lower.limits = self$lower.limits,
        upper.limits = self$upper.limits,
        # maxit = self$maxit,
        type.gaussian = self$type.gaussian,
        type.logistic = self$type.logistic,
        standardize.response = self$standardize.response,
        type.multinomial = self$type.multinomial,
        relax = self$relax,
        # trace.it = self$trace.it,
        cox.ties = self$cox.ties,
        control = self$control,
        s = self$s
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
