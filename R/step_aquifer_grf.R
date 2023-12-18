#' R6 Class
#'
#' `StepAquiferGRF` Generates the drawdown using the Generalized Radial Flow (GRF)
#' model. This method defaults to a fast FFT convolution so many rates can be included,
#' but requires a regular time series.
#'
#' @param time numeric vector of input times
#' @param flow_rate numeric vector of flow rates for each input time
#' @param thickness numeric thickness of aquifer
#' @param radius numeric radial distance to monitoring well (> 0.0)
#' @param specific_storage the specific storage
#' @param hydraulic_conductivity the hydraulic conductivity
#' @param flow_dimension numeric value for the flow dimension: 1 = linear,
#' 2 = radial (i.e. Theis), 3 = spherical.
#'
#' @inheritParams Step
#'
#' @family aquifer
#'
#' @export
StepAquiferGRF <- R6Class(

  classname = 'step_aquifer_grf',
  inherit = Step,

  public = list(
    # step specific variables
    time = NULL,
    flow_rate = NULL,
    thickness = NULL,
    radius = NULL,
    specific_storage = NULL,
    hydraulic_conductivity = NULL,
    flow_dimension = NULL,

    initialize = function(time,
                          flow_rate,
                          thickness = 1.0,
                          radius = 100.0,
                          specific_storage = 1.0e-6,
                          hydraulic_conductivity = 1.0e-4,
                          flow_dimension = 2.0,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {

      # get function parameters to pass to parent
      step_name    <- "step_aquifer_grf"
      type         <- 'add'
      enq <- rlang::enquos(time, flow_rate)
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize))
      )
      do.call(super$initialize, inputs)

      # step specific values
      self$time = enquos(time)
      self$flow_rate = enquos(flow_rate)
      self$thickness = thickness
      self$radius = radius
      self$specific_storage = specific_storage
      self$hydraulic_conductivity = hydraulic_conductivity
      self$flow_dimension = flow_dimension

      invisible(self)
    },

    bake = function(new_data) {
      setNames(grf_time(radius = self$radius,
                        specific_storage = self$specific_storage,
                        hydraulic_conductivity = self$hydraulic_conductivity,
                        thickness = self$thickness,
                        time = new_data[[1]],
                        flow_rate = new_data[[2]],
                        flow_dimension = self$flow_dimension), self$id)
    }

  )
)



