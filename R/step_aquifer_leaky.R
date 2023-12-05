#' R6 Class
#' @title
#' hantush_jacob
#'
#' @description
#' Convolution of hantush well function and flow rates in the frequency domain.
#' Time series needs to be regularily spaced.
#'
#' @param radius distance to monitoring interval
#' @param storativity aquifer storativity
#' @param transmissivity aquifer transmissivity
#' @param leakage hantush leakage
#' @param time prediction times
#' @param flow_rate well flow rates
#' @param n_terms number of terms to use in Hantush solution.  More is more precise but slower.
#' @inheritParams Step
#'
#' @family aquifer
#'
#' @export
StepAquiferLeaky <- R6Class(

  classname = 'step_aquifer_leaky',
  inherit = Step,

  public = list(

    # step specific variables
    time = NULL,
    flow_rate = NULL,
    leakage = NULL,
    radius = NULL,
    storativity = NULL,
    transmissivity = NULL,
    max_terms = NULL,

    initialize = function(time,
                          flow_rate,
                          leakage = 100,
                          radius = 100,
                          storativity = 1e-5,
                          transmissivity = 1e-4,
                          max_terms = 20,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {

      # get function parameters to pass to parent
      step_name    <- "step_aquifer_leaky"
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
      self$flow_rate = enquos(flow_rate)

      self$leakage = leakage
      self$radius = radius

      # K & Ss
      self$storativity = storativity
      self$transmissivity = transmissivity
      self$max_terms = max_terms

      invisible(self)
    },

    bake = function(new_data) {
      setNames(
        list(hantush_jacob(
          new_data[[1]],
          new_data[[2]],
          self$radius,
          self$storativity,
          self$transmissivity,
          self$leakage,
          self$max_terms
        )),
        self$id)
    }
  )
)


