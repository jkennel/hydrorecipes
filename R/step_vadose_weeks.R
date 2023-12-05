#' R6 Class
#'
#' `StepVadoseWeeks` Weeks solution
#'
#' @param time
#' @param flow_rate
#' @param thickness
#' @param radius
#' @param radius_patch
#' @param hydraulic_conductivity_inner
#' @param hydraulic_conductivity_outer
#' @param specific_storage_inner
#' @param specific_storage_outer
#'
#' @inheritParams Step
#'
#' @family vadose
#'
#' @export
StepVadoseWeeks <- R6Class(

  classname = 'step_vadose_weeks',
  inherit = Step,

  public = list(

    # step specific variables
    time = NULL,
    air_diffusivity = NULL,
    thickness = NULL,
    precision = NULL,
    inverse = NULL,

    initialize = function(time,
                          air_diffusivity = 0.2,
                          thickness = 40.0,
                          precision = 1e-12,
                          inverse = FALSE,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {

      # get function parameters to pass to parent
      step_name    <- "step_vadose_weeks"
      type         <- 'add'
      enq <- rlang::enquos(time)
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize))
      )
      do.call(super$initialize, inputs)

      # step specific values
      self$time = enquos(time)
      self$air_diffusivity = air_diffusivity
      self$thickness = thickness
      self$precision = precision
      self$inverse = inverse

      invisible(self)
    },

    bake = function(new_data) {
      setNames(
        list(vadose_response(
          new_data[[1]],
          self$air_diffusivity,
          self$thickness,
          self$precision,
          self$inverse,
        )),
        self$id)
    }
  )
)



