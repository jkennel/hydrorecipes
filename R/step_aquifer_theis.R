#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Theis Step -------------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepAquiferTheis` Generates the drawdown using the Generalized
#' Radial Flow (GRF) model with flow dimension equal to 2. This method defaults
#' to a fast FFT convolution so many rates can be included, but requires a
#' regular time series.
#'
#' @inheritParams StepAquiferGRF
#'
#' @family aquifer
#'
#' @references
#' Barker, J.A., A generalized radial flow model for hydraulic tests
#'  in fractured rock. Water Resour. Res., 24 (1988), pp. 1796-1804,
#'  10.1029/WR024i010p01796
#'
#' Theis, C.V., 1935: The relation between the lowering of the piezometric
#'  surface and the rate and duration of discharge of a well using
#'  ground-water storage, Transactions of the American Geophysical Union,
#'  16th Annual Meeting, Part 2, pp. 519-524.
#'
#' @export
StepAquiferTheis <- R6Class(
  classname = "step_aquifer_theis",
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
                          ...) {
      # get function parameters to pass to parent
      inputs <- list(
        time = substitute(time),
        flow_rate = substitute(flow_rate),
        thickness = thickness,
        radius = radius,
        specific_storage = specific_storage,
        hydraulic_conductivity = hydraulic_conductivity,
        role = role,
        ...
      )

      do.call(super$initialize, inputs)

      invisible(self)
    }
  )
)
