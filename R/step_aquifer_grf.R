#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Generalized Radial Flow (GRF) Step -------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepAquiferGRF` Generates the drawdown using the Generalized
#' Radial Flow (GRF) model. This method defaults to a fast FFT
#' convolution so many rates can be included, but requires a regular
#' time series.
#'
#' @inheritParams Step
#'
#' @references
#' Barker, J.A., A generalized radial flow model for hydraulic tests
#'  in fractured rock. Water Resour. Res., 24 (1988), pp. 1796-1804,
#'  10.1029/WR024i010p01796
#'
#' @family aquifer
#'
#' @export
StepAquiferGRF <- R6Class(
  classname = "step_aquifer_grf",
  inherit = Step,
  public = list(
    # step specific variables
    #' @field time numeric vector of input times
    time = NULL,
    #' @field flow_rate numeric vector of flow rates for each input time
    flow_rate = NULL,
    #' @field thickness numeric thickness of aquifer
    thickness = NULL,
    #' @field radius numeric radial distance to monitoring well (> 0.0)
    radius = NULL,
    #' @field specific_storage the specific storage
    specific_storage = NULL,
    #' @field hydraulic_conductivity the hydraulic conductivity
    hydraulic_conductivity = NULL,
    #' @field flow_dimension numeric value for the flow dimension: 1 = linear,
    #' 2 = radial (i.e. Theis), 3 = spherical.
    flow_dimension = NULL,
    initialize = function(time,
                          flow_rate,
                          thickness = 1.0,
                          radius = 100.0,
                          specific_storage = 1.0e-6,
                          hydraulic_conductivity = 1.0e-4,
                          flow_dimension = 2.0,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      flow_rate <- deparse(substitute(flow_rate))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_aquifer_grf"
      env_list$type <- "add"
      super$initialize(
        terms = c(as.symbol(time), as.symbol(flow_rate)),
        env_list
      )


      # step specific values
      self$time <- time
      self$flow_rate <- flow_rate
      self$thickness <- thickness
      self$radius <- radius
      self$specific_storage <- specific_storage
      self$hydraulic_conductivity <- hydraulic_conductivity
      self$flow_dimension <- flow_dimension

      self$columns <- c(time, flow_rate)

      invisible(self)
    },
    bake = function(new_data) {
      setNames(grf_time(
        radius = self$radius,
        specific_storage = self$specific_storage,
        hydraulic_conductivity = self$hydraulic_conductivity,
        thickness = self$thickness,
        time = new_data[[self$time]],
        flow_rate = new_data[[self$flow_rate]],
        flow_dimension = self$flow_dimension
      ), self$id)
    }
  )
)
