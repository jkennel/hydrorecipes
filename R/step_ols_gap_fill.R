#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Fill in Gaps using Regression ------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepOlsGapFill`
#'
#' @param terms
#' @param recipe Recipe to use for filling gaps
#'
#' @inheritParams Step
#'
#' @family gap_fill
#'
#' @export
StepOlsGapFill <- R6Class(
  classname = "step_ols_gap_fill",
  inherit = Step,
  public = list(

    # step specific variables
    recipe = NULL,
    coefficients = NULL,
    initialize = function(terms,
                          recipe,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_ols_gap_fill"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )


      # step specific values
      self$recipe <- recipe

      invisible(self)
    },
    bake = function(new_data) {

      rec <- self$recipe
      rec <- rec$prep()$bake(data = new_data)
      dat <- rec$plate(type = "list")

      x <- get_regression_data(dat, rec$term_info, id_type = "predictor")
      y <- get_regression_data(dat, rec$term_info, id_type = "outcome")

      mode(x$data) <- "double"
      mode(y$data) <- "double"

      # ti <- collapse::qDF(rec$term_info)
      # ti <- ti[ti$source != "removed", ]



      # outcomes <- ti[ti$roles == "outcome", ]
      # predictors <- ti[ti$roles == "predictor", ]
      #
      # outcome_ids <- which(nms %in% outcomes$variable)
      # predictor_ids <- which(nms %in% predictors$variable)
      #
      # m_predictors <- collapse::qM(unclass(dat)[predictor_ids])
      # m_outcomes <- collapse::qM(unclass(dat)[outcome_ids])

      # remove na values in the outcomes
      self$coefficients <- determine_coefficients(x, y)


      lst <- collapse::mctl(x$data[, , drop = FALSE] %*% self$coefficients[, , drop = FALSE])

      self$new_columns <- name_columns(self$prefix, colnames(y$data), ncol(y$data))
      names(lst) <- self$new_columns

      lst

      # rec <- self$recipe
      # rec <- rec$prep()$bake(data = new_data)
      # ti <- collapse::qDF(rec$term_info)
      # ti <- ti[ti$source != "removed", ]
      # dat <- rec$plate(type = "list")
      # nms <- names(dat)
      #
      #
      #
      # outcomes <- ti[ti$roles == "outcome", ]
      # predictors <- ti[ti$roles == "predictor", ]
      #
      # outcome_ids <- which(nms %in% outcomes$variable)
      # predictor_ids <- which(nms %in% predictors$variable)
      #
      # m_predictors <- collapse::qM(unclass(dat)[predictor_ids])
      # m_outcomes <- collapse::qM(unclass(dat)[outcome_ids])
      # mode(m_predictors) <- "double"
      # mode(m_outcomes) <- "double"
      #
      # # remove na values in the outcomes
      # wh <- which(!is.na(m_outcomes))
      #
      # # solve
      # fit <- llt_solve(
      #   m_predictors[wh, , drop = FALSE],
      #   m_outcomes[wh, , drop = FALSE]
      # )
      # self$coefficients <- fit
      #
      # lst <- collapse::mctl(m_predictors[, , drop = FALSE] %*% fit[, , drop = FALSE])
      #
      # self$new_columns <- name_columns(self$prefix, outcomes$variable, length(outcome_ids))
      # names(lst) <- self$new_columns
      #
      # lst
    }
  )
)




# set.seed(123)
# frm <- formula(x ~ y + z)
# x <- cumsum(rnorm(100))
# dat <- data.table(x = x, y = x, z = as.numeric(1:100))
# dat[, x := x + c(rep(20, 50), rep(0, 50))]
# dat[, x := x + 3 *sin(z * 0.1)]
# tmp <- copy(dat$x)
# plot(tmp, type = 'l')
# dat[60:70, x := NA_real_]
#
#
# dat <- unclass(dat)
# f = Recipe$new(formula = frm, data = (dat))$
#   add_step(StepFindInterval$new(z, vec = c(0, 50.5, 101)))$
#   add_step(StepIntercept$new())$
#   add_step(StepSplineB$new(z, df = 10))$
#   # add_step(StepDropColumns$new(y))$
#   add_step(StepDropColumns$new(z))
#
# frec = Recipe$new(formula = frm, data = dat)$
#       add_step(StepRegressionGapFill$new(c(x, y, z), recipe = f))
#
# tmp <- frec$prep()$bake()$plate()
# points(dat$x, type = 'p', pch = 20)
# points(tmp$update, type = 'l', col = 'red')
