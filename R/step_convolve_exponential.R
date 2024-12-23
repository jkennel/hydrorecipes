#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# FFT Convolution of a Term with a Gamma Kernel --------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepConvolveExponential <- R6Class(
  classname = "step_convolve_exponential",
  inherit = StepConvolveGamma,
  public = list(

    initialize = function(terms,
                          amplitude,
                          theta,
                          align = "right",
                          max_length = Inf,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      inputs <- list(
        terms = substitute(terms),
        amplitude = amplitude,
        k = 1.0,
        theta = theta,
        align = align,
        max_length = max_length,
        role = role,
        step_name = "step_convolve_exponential",
        ...
      )

      do.call(super$initialize, inputs)

      invisible(self)
    }
  )
)
