#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# FFT Convolution of a Term with a Gamma Kernel --------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepConvolveGamma <- R6Class(
  classname = "step_convolve_gamma",
  inherit = StepKernelFilter,
  public = list(

    amplitude = NULL,
    k = NULL,
    theta = NULL,
    cutoff = NULL,

    initialize = function(terms,
                          amplitude,
                          k,
                          theta,
                          align = "right",
                          max_length = Inf,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      inputs <- list(
        terms = substitute(terms),
        kernel = NA_real_,
        align = align,
        role = role,
        step_name = "step_convolve_gamma",
        ...
      )

      do.call(super$initialize, inputs)

      self$amplitude <- amplitude
      self$k <- k
      self$theta <- theta

      if (is.finite(max_length)) {
        self$n_na_max <- max_length
      } else {
        if (!is.null(self$varying)) {
          if (!is.null(self$varying$upper)) {
            self$n_na_max <- 720L * ceiling(self$varying$upper[3])
          }
        } else {
          self$n_na_max <- 720L * ceiling(self$theta)
        }
      }


      invisible(self)
    },

    bake = function(s) {

      self$kernel <- list(gamma_3(0L:(self$n_na_max),
                                  self$amplitude,
                                  self$k,
                                  self$theta))

      super$bake(s)

    }
  )
)
