#' R6 Class
#'
#' `StepBSpline` generates basis splines.
#'
#' @param df
#' @param internal_knots
#' @param degree
#' @param intercept
#' @param boundary_knots
#' @param periodic
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
                          degree = 3L,
                          intercept = FALSE,
                          boundary_knots = NULL,
                          periodic = FALSE,
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
    bake = function(new_data) {

      column_name     <- self$columns

      basis <- b_spline_list(x = unclass(new_data)[[column_name]],
                         df = self$df,
                         degree = 3L,
                         internal_knots = self$internal_knots,
                         boundary_knots = self$boundary_knots
                         )

      names(basis) <- file.path(self$id, seq_len(length(basis)), column_name, fsep = '_')
      self$result <- append(self$result, basis)

      self$result
    }
  )
)



