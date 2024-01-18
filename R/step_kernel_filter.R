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

    initialize = function(terms,
                          kernel,
                          align = "center",
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- 'step_kernel_filter'
      env_list$type <- 'add'
      super$initialize(terms = terms,
                       env_list[names(env_list) != "terms"])


      # step specific values
      self$kernel <- if (!inherits(kernel, "list")) list(kernel) else kernel
      n_kernel <- length(kernel)
      n_align <- length(align)
      if (n_align != 1) {
        stop('align should be length 1')
      }
      self$align  <- align

      invisible(self)
    },

    bake = function(new_data) {

      column_name     <- self$columns

      filt <- list()
      for (i in seq_along(column_name)) {

          if (self$align == "center") {
            filt[[i]] <- convolve_overlap_save_list(unclass(new_data)[[i]],
                                                    self$kernel, 1)
          }
          if (self$align == "right") {
            filt[[i]] <- convolve_overlap_save_list(unclass(new_data)[[i]],
                                                    self$kernel, 0)
          }
          if (self$align == "left") {
            filt[[i]] <- convolve_overlap_save_list(unclass(new_data)[[i]],
                                                    self$kernel, 2)
          }

        names(filt[[i]]) <- name_columns(self$id, column_name, length(self$kernel))
      }

      unlist(filt, recursive = FALSE)

    }

  )
)



