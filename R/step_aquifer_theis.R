#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Theis Step -------------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
        step_name = "step_aquifer_theis",
        ...
      )

      do.call(super$initialize, inputs)

      invisible(self)
    }
  )
)


