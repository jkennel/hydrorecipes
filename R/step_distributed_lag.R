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
    # lag_max integer the maximum lag.
    lag_max = NULL,
    # intercept boolean use an intercept for the basis matrix
    intercept = NULL,
    # basis_matrix matrix the basis matrix.
    basis_matrix = NULL,


    initialize = function(terms,
                          n_lag = 12L,
                          lag_max = 86400L,
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
          self$knots <- log_lags(self$n_lag, self$lag_max)
        }

        self$n_lag <- length(knots)
        self$lag_max <- max(knots)
        self$intercept <- intercept

        rng = 0:self$lag_max
        one_n = c(1L, self$n_lag)

        # natural spline
        self$basis_matrix <- n_spline_list(rng, 0L, 3L, self$knots[-one_n],
                                           self$knots[one_n], self$intercept,
                                           FALSE, 0L, FALSE)

      } else {
        self$lag_max <- nrow(basis_matrix)
        self$n_lag <- ncol(basis_matrix)
        self$basis_matrix <- collapse::mctl(basis_matrix)
      }

      self$n_na_max <- self$lag_max

      invisible(self)
    },
    bake = function(s) {

      # column_name <- self$columns
      dl <- list()
      for (i in seq_along(self$columns)) {
        column_name <- self$columns[i]
        dl[[i]] <- distributed_lag_list4(
          s[["result"]][[column_name]],
          self$basis_matrix,
          self$lag_max
        )

        names(dl[[i]]) <- name_columns(self$prefix, column_name[i], length(dl[[i]]))
        self$new_columns <- c(self$new_columns, names(dl[[i]]))
      }
      names(self$basis_matrix) <- self$new_columns
      # self$columns <- rep(self$columns, each = length(self$basis_matrix))

      self$result <- unlist(dl, recursive = FALSE)
      # self$result
      return(NULL)
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
      resp <- basis_matrix %*% co

      # print("mult")
      # print(str(basis_matrix[, wh, drop = FALSE]))
      # print(str(co[wh, , drop = FALSE]))
      list(x = rep(0:(nr - 1L), nc * 2L),
           variable = rep(c("coefficient", "cumulative"), each = nr * nc),
           value = c(resp, collapse::fcumsum(resp)),
           step_id = rep.int(self$id, 2L * nr * nc),
           outcome = rep(rep(colnames(co), each = nr), 2L),
           term = rep.int(self$prefix, 2L * nr * nc),
           step_columns = rep.int(paste(self$columns, collapse = "_"), 2L * nr * nc))
    }
  )
)
