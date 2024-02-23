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
    #' @field n_lag integer the maximum lag.
    max_lag = NULL,
    initialize = function(terms,
                          knots,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_distributed_lag"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"]
      )

      # step specific values
      self$knots <- knots
      self$n_lag <- length(knots)
      self$max_lag <- max(knots)
      invisible(self)
    },
    bake = function(new_data) {
      column_name <- self$columns

      dl <- list()
      for (i in seq_along(column_name)) {
        dl[[i]] <- distributed_lag_list3(
          unclass(new_data)[[i]],
          self$n_lag,
          self$max_lag + 1,
          0L,
          3L,
          self$knots[2:(self$n_lag - 1L)],
          self$knots[c(1, self$n_lag)],
          TRUE,
          FALSE,
          0L,
          FALSE
        )
        names(dl[[i]]) <- name_columns(self$id, column_name[i], length(dl[[i]]))
      }
      unlist(dl, recursive = FALSE)
    }
  )
)
