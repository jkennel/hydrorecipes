#' R6 Class
#'
#' `StepLeadLag` generates lagged (or leading) vectors.
#'
#' @param lag integer vector for the amount to lag or lead.  A negative value indicates leading.
#' @param n_shift the amount to shift the starting point. When `n_subset` is not
#' equal to 1 it may be useful to start the lagging at an offset from the beginning
#' of the series. Default is 0.
#' @param n_subset the spacing in rows between values.
#'
#' @inheritParams Step
#'
#' @export
StepLeadLag <- R6Class(
  classname = 'step_lead_lag',
  inherit = Step,

  public = list(

    # step specific variables
    lag = NULL,
    n_shift = NULL,
    n_subset = NULL,

    initialize = function(...,
                          lag,
                          n_shift = 0L,
                          n_subset = 1L,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_lead_lag"
      type         <- 'add'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)

      # step specific values
      self$lag      <- as.integer(sort(lag))
      self$n_shift  <- as.integer(n_shift)
      self$n_subset <- as.integer(n_subset)

      invisible(self)
    },

    bake = function(new_data) {

      column_name <- self$columns

      if(self$n_subset == 1) {
        ll <- collapse::flag(list(new_data), self$lag)
      } else {
        ll <- lag_list(new_data,
                       self$lag,
                       n_subset = self$n_subset,
                       n_shift = self$n_shift)
      }

      names(ll) <- file.path(self$id,
                             column_name,
                             pad_num(length(self$lag)),
                             fsep = '_')

      ll

    }

  )
)



