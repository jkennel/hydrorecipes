#' R6 Class
#'
#' `StepKernelFilter` linearly convolve a kernel with a data series.
#'
#' @param kernel the convolution kernel
#'
#' @inheritParams Step
#'
#' @export
StepKernelFilter <- R6Class(
  classname = 'step_kernel_filter',
  inherit = Step,

  public = list(

    # step specific variables
    kernel = NULL,
    align = NULL,

    initialize = function(...,
                          kernel,
                          align = "center",
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name <- "step_kernel_filter"
      type      <- 'add'
      enq       <- NULL
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)

      # step specific values
      self$kernel <- kernel
      self$align  <- align

      invisible(self)
    },

    bake = function(new_data) {

      if(self$align == "center") {
        new_data <- convolve_overlap_save_list(new_data, self$kernel, 1)
      }
      if(self$align == "right") {
        new_data <- convolve_overlap_save_list(new_data, self$kernel, 0)
      }
      if(self$align == "left") {
        new_data <- convolve_overlap_save_list(new_data, self$kernel, 2)
      }

      names(new_data) <- file.path(self$id,
                                   pad_num(length(new_data)),
                                   fsep = "_")

      return(new_data)

    }

  )
)



