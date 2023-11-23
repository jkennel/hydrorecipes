#' R6 Class
#'
#' `StepBSpline` generates basis splines.
#'
#' @param df
#' @param internal_knots
#' @param boundary_knots
#' @param intercept
#' @param periodic
#' @param degree
#' @inheritParams Step
#'
#' @export
StepBSpline <- R6Class(
  classname = 'step_b_spline',
  inherit = Step,

  public = list(

    # step specific variables
    df = 0L,
    internal_knots = NULL,
    degree = 3L,
    intercept = FALSE,
    boundary_knots = NULL,
    periodic = FALSE,

    initialize = function(...,
                          df = 0L,
                          internal_knots = NULL,
                          boundary_knots = NULL,
                          intercept = FALSE,
                          periodic = FALSE,
                          degree = 3L,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_b_spline"
      type         <- 'add'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1])
      )
      do.call(super$initialize, inputs)

      # step specific values
      self$df <- df
      self$internal_knots  <- as.numeric(sort(internal_knots))
      self$degree <- degree
      self$intercept <- intercept
      self$boundary_knots <- as.numeric(boundary_knots)
      self$periodic <- periodic

      invisible(self)

    },
    prep = function(new_data) {

      if (self$df != 0L) {
        ik <- collapse::fquantile(new_data,
                                  probs = seq(0, 1, self$df),
                                  na.rm = TRUE)
        self$boundary_knots = ik[c(1, length(ik))]
        self$internal_knots <- ik[-c(1, length(ik))]
      }

    },
    bake = function(new_data) {

      column_name     <- self$columns

      basis <- b_spline_list(x = new_data,
                             df = self$df,
                             degree = 3L,
                             internal_knots = self$internal_knots,
                             boundary_knots = self$boundary_knots
      )

      names(basis) <- file.path(self$id,
                                seq_len(length(basis)),
                                column_name,
                                fsep = '_')

      basis
    }
  )
)



