#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Distributed Lag Step ---------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepDistributedLag <- R6Class(

  classname = "step_distributed_lag",
  inherit = Step,
  public = list(

    # step specific variables
    # knots the locations of the knots for the basis matrix.
    knots = NULL,
    # n_lag integer the number of lag terms.
    n_lag = NULL,
    # max_lag integer the maximum lag.
    max_lag = NULL,
    # intercept boolean use an intercept for the basis matrix
    intercept = NULL,
    # basis_matrix matrix the basis matrix.
    basis_matrix = NULL,


    initialize = function(terms,
                          n_lag = 12L,
                          max_lag = 86400L,
                          knots = NA_real_,
                          basis_matrix = NA_real_,
                          intercept = FALSE,
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
        } else {
          self$knots <- log_lags(self$n_lag, self$max_lag)
        }

        self$n_lag <- length(knots)
        self$max_lag <- max(knots)
        self$intercept <- intercept

        rng = 0:self$max_lag
        one_n = c(1L, self$n_lag)

        self$basis_matrix <- n_spline_list(rng, 0L, 3L, self$knots[-one_n],
                                           self$knots[one_n], self$intercept,
                                           FALSE, 0L, FALSE)

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
      names(self$basis_matrix) <- self$new_columns
      # self$columns <- rep(self$columns, each = length(self$basis_matrix))

      unlist(dl, recursive = FALSE)

    },
    # returns a named list
    response = function(co) {

      basis_matrix <- collapse::qM(self$basis_matrix)


      nr <- nrow(basis_matrix)
      nc <- ncol(co)

      if (nrow(co) != ncol(basis_matrix)) {
        warning("The provided formula does not match the basis_matrix.
              Did you select fewer distributed lag terms?")
      }

      wh <- intersect(colnames(basis_matrix), rownames(co))

      # check for multiple outcomes!!
      resp <- basis_matrix[, wh, drop = FALSE] %*% co[wh, , drop = FALSE]

      list(x = rep(0:(nr - 1L), nc * 2L),
           variable = rep(c("coefficient", "cumulative"), each = nr * nc),
           value = c(resp, collapse::fcumsum(resp)),
           step_id = rep(self$id, 2L * nr * nc),
           outcome = rep(rep(colnames(co), each = nr), 2L),
           term = rep("distributed_lag_interpolated", 2L * nr * nc))
    }
  )
)
