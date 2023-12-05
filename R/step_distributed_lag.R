#' R6 Class
#'
#' `StepDistributedLag` generates distributed lag vectors.
#'
#' @param knots
#' @inheritParams Step
#'
#' @export
StepDistributedLag <- R6Class(
  classname = 'step_distributed_lag',
  inherit = Step,

  public = list(

    # step specific variables
    knots = NULL,
    n_lag = NULL,
    max_lag = NULL,

    initialize = function(...,
                          knots,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_distributed_lag"
      type         <- 'add'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)

      # step specific values
      self$knots <- knots
      self$n_lag <- length(knots)
      self$max_lag <- max(knots)
      invisible(self)
    },

    bake = function(new_data) {

      column_name <- self$columns

      dl <- distributed_lag_list3(
        new_data,
        self$n_lag,
        self$max_lag,
        0L,
        3L,
        self$knots[2:(self$n_lag - 1L)],
        self$knots[c(1, self$n_lag)],
        TRUE,
        FALSE,
        0L,
        FALSE
      )

      names(dl) <- file.path(self$id,
                             column_name,
                             pad_num(length(dl)),
                             fsep = '_')

      dl

    }

  )
)
