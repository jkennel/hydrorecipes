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
    max_length = NULL,

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
      self$max_length <- max_length

      invisible(self)
    },

    bake = function(new_data) {

      n <- length(unclass(new_data)[[1]])

      max_t <- 720L * ceiling(self$theta);

      self$max_length <- pmin(n - 1L, self$max_length)
      max_t <- pmin(max_t, self$max_length)

      self$kernel <- list(gamma_3(1L:max_t, self$amplitude, self$k, self$theta))

      super$bake(new_data)

    }
  )
)
