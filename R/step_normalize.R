#' R6 Class
#'
#' `StepNormalize` adjust the dispersion by standard deviation.
#' @inheritParams Step
#'
#' @export
StepNormalize <- R6Class(
  classname = 'step_normalize',
  inherit = Step,

  public = list(
    center = c(),
    scale = c(),
    na_rm = NA,
    fun = NULL,
    # step specific variables
    initialize = function(...,
                          role = "predictor",
                          skip = FALSE,
                          na_rm = TRUE,
                          keep_original_cols = FALSE) {

      # get function parameters to pass to parent
      step_name    <- "step_normalize"
      type         <- 'modify'
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize)[-1L])
      )
      do.call(super$initialize, inputs)

      self$na_rm <- na_rm

      invisible(self)
    },
    prep = function(new_data) {

      self$center <- collapse::fmean(new_data,
                                     na.rm = self$na_rm)
      self$scale <- collapse::fsd(new_data,
                                  na.rm = self$na_rm)
    },
    # subtract the central value from a column
    bake = function(new_data) {

      for(i in seq_along(self$columns)) {
        new_data[[i]] %-=% self$center[i]
        new_data[[i]] %*=% (1.0 / self$scale[i])
      }

      # fscale(new_data, self$center, self$scale)
      # scale_list_param_std(new_data,
      #                        center = self$center,
      #                        scale = 1.0/self$scale)

      new_data
    }
  )
)

