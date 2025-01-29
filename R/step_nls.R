#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Predict Regression Terms -----------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepNls <- R6Class(
  classname = "step_nls",
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
    # start = NULL,
    # lower = NULL,
    # upper = NULL,
    algorithm = NULL,
    n_subset = NULL,
    n_shift = NULL,
    range = NULL,
    control = NULL,

    trace = NULL,
    # do_predict = NULL,
    # s = NULL,
    # df_residual = NULL,
    # rank = NULL,
    # std_error = NULL,

    initialize = function(formula = NULL,
                          # start = 0,
                          # lower = 0,
                          # upper = 0,
                          algorithm = "lm",
                          n_subset = 1L,
                          n_shift = 0L,
                          range = c(-Inf, Inf),
                          control =  gsl_nls_control(xtol = 1e-8),
                          trace = FALSE,
                          role = "predictor",
                          # do_response = TRUE,
                          # do_predict = TRUE,
                          ...) {

      require(gslnls)
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_nls"
      env_list$type <- "model"
      super$initialize(
        terms = NULL,
        env_list[names(env_list) != "terms"],
        ...
      )

      self$formula   <- formula
      self$algorithm <- algorithm
      self$n_subset  <- n_subset
      self$n_shift   <- n_shift
      self$range     <- range
      self$control   <- control
      self$trace     <- trace
      # self$do_response <- do_response
      # self$do_predict <- do_predict

      invisible(self)
    },
    fn = function(pars, r, varying, inds) {

      # update all steps
      r$update_steps(pars, varying)
      nc <- r$update_pars(pars, varying)

      # get step numbers with varying parameters
      s  <- unique(varying[["step_number"]])

      # get step numbers with coefficient parameters
      step_dont_rerun <- which(r$get_step_is_coef())

      # multiply list columns by coefficients and subset
      list_multiply_subset(r$plate(type = "list",
                                   steps = s,
                                   step_dont_rerun = step_dont_rerun), nc, inds - 1L)

    },
    bake = function(r) {

      # get parameters to vary
      varying <- r$get_varying()

      y      <- r$get_outcome_variable()[[1L]]
      n_y    <- length(y)
      to_rem <- max(unlist(r$get_step_field("n_na_max")), na.rm = TRUE) + 1L

      # set range for residuals
      if (is.null(self$range)) {
        inds   <- seq(to_rem + self$n_shift, n_y, self$n_subset)
      } else {
        self$range[1L] <- pmax(to_rem, self$range[1L])
        self$range[2L] <- pmin(self$range[2L], n_y)
        inds          <- seq(self$range[1L] + self$n_shift,
                             self$range[2L],
                             self$n_subset)
      }

      # print(collapse::qDF(varying))
      # print(self$range)

      self$fit <- gsl_nls(

        fn = self$fn,
        y = y[inds],

        start = c(pars = varying[["start"]]),
        lower = c(pars = varying[["lower"]]),
        upper = c(pars = varying[["upper"]]),

        control = self$control,
        trace = self$trace,

        algorithm = self$algorithm,
        r = r,
        varying = varying,
        inds = inds
      )

      # print(self$fit)


      return(NULL)


    }
  )
)
