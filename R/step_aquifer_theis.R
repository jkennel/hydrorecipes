#' R6 Class
#'
#' `StepAquiferTheis` generates lagged (or leading) vectors.
#'
#' @inheritParams StepAquiferGRF
#'
#' @family aquifer
#'
#' @export
StepAquiferTheis <- R6Class(

  classname = 'step_aquifer_theis',
  inherit = StepAquiferGRF,

  public = list(
    # step specific variables

    initialize = function(time,
                          flow_rate,
                          thickness = 1.0,
                          radius = 100.0,
                          specific_storage = 1.0e-6,
                          hydraulic_conductivity = 1.0e-4,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {

      step_name    <- "step_aquifer_theis"
      type         <- 'add'
      flow_dimension = 2.0

      time = enquo(time)
      flow_rate = enquo(flow_rate)

      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize))
      )
      do.call(super$initialize, inputs)

      invisible(self)
    }

  )
)



