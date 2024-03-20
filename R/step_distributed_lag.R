#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Distributed Lag Step ---------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepDistributedLag` generates distributed lag vectors.
#'
#'
#' @inheritParams Step
#'
#' @export
StepDistributedLag <- R6Class(
  classname = "step_distributed_lag",
  inherit = Step,
  public = list(

    # step specific variables
    #' @field knots the locations of the knots for the basis matrix.
    knots = NULL,
    #' @field n_lag integer the number of lag terms.
    n_lag = NULL,
    #' @field max_lag integer the maximum lag.
    max_lag = NULL,
    #' @field basis_matrix matrix the basis matrix.
    basis_matrix = NULL,


    initialize = function(terms,
                          n_lag = 12L,
                          max_lag = 86400L,
                          knots = NA_real_,
                          basis_matrix = NA_real_,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_distributed_lag"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"],
        ...
      )

      # set up basis matrix
      if (all(is.na(basis_matrix))) {
        # step specific values
        if (!all(is.na(knots))) {
          self$knots <- knots
          self$n_lag <- length(knots)
          self$max_lag <- max(knots)
        } else {
          self$knots <- log_lags_arma(self$n_lag, self$max_lag)
          self$n_lag <- length(knots)
          self$max_lag <- max(knots)
        }

        rng = 0:self$max_lag
        one_n = c(1L, self$n_lag)

        self$basis_matrix <- n_spline_list(rng, 0L, 3L, self$knots[-one_n],
                                           self$knots[one_n], FALSE, FALSE,
                                           0L, FALSE)

      } else {
        self$max_lag <- nrow(basis_matrix)
        self$n_lag <- ncol(basis_matrix)
        self$basis_matrix <- collapse::mctl(basis_matrix)
      }


      invisible(self)
    },
    bake = function(new_data) {

      column_name <- self$columns
      dl <- list()
      for (i in seq_along(column_name)) {

        dl[[i]] <- distributed_lag_list4(
          unclass(new_data)[[i]],
          self$basis_matrix,
          self$max_lag
        )

        names(dl[[i]]) <- name_columns(self$prefix, column_name[i], length(dl[[i]]))
        self$new_columns <- c(self$new_columns, names(dl[[i]]))
      }

      self$columns <- rep(self$columns, each = length(self$basis_matrix))

      unlist(dl, recursive = FALSE)

    },
    response = function(co) {

      basis_matrix <- collapse::qM(self$basis_matrix)
      n <- nrow(basis_matrix)

      # print(str(basis_matrix))
      # print(co)
      # check for multiple outcomes!!
      resp <- basis_matrix %*% co

      variable <- c(
        rep("coefficient", n),
        rep("cumulative", n)
      )

      value <- c(
        resp,
        cumsum(resp)
      )

      list(x = rep(0:(n-1), 2L),
           variable = variable,
           value = value,
           step_id = rep(self$id, 2L * n),
           term = "distributed_lag_interpolated")
    }
  )
)
