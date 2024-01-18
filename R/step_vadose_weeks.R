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
                          ...) {

      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      env_list <- get_function_arguments()
      env_list$step_name <- 'step_vadose_weeks'
      env_list$type <- 'add'
      super$initialize(terms = c(as.symbol(time)),
                       env_list)

      # step specific values
      self$time = time
      self$air_diffusivity = air_diffusivity
      self$thickness = thickness
      self$precision = precision
      self$inverse = inverse
      self$columns <- time

      invisible(self)
    },

    bake = function(new_data) {

        vadose_response(
          new_data[[1]],
          self$air_diffusivity,
          self$thickness,
          self$precision,
          self$inverse
        )
    }
  )
)



