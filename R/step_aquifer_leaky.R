#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Hantush and Jacob Leaky Aquifer Step -----------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#' @title
#' hantush_jacob 1955 solution for a leaky aquifer
#'
#' @description
#' Convolution of hantush well function and flow rates in the frequency domain.
#' Time series needs to be regularly spaced.
#'
#'
#' @inheritParams Step
#'
#' @references
#' J.H.A. Prodanoff; W.J. Mansur; F.C.B. Mascarenhas (2006).
#'  Numerical evaluation of Theis and Hantush-Jacob well functions. , 318(1-4),
#'  0–183. doi:10.1016/j.jhydrol.2005.05.026 eq: 10, 11, 12
#'
#' Hantush, M.S. and C.E. Jacob, 1955. Non-steady radial flow in an infinite
#'  leaky aquifer, Am. Geophys. Union Trans., vol. 36, no. 1, pp. 95-100.
#'
#' @family aquifer
#'
#' @export
StepAquiferLeaky <- R6Class(
  classname = "step_aquifer_leaky",
  inherit = Step,
  public = list(

    # step specific variables
    #' @field time prediction times
    time = NULL,
    #' @field flow_rate well flow rates
    flow_rate = NULL,
    #' @field leakage hantush leakage
    leakage = NULL,
    #' @field radius distance to monitoring interval
    radius = NULL,
    #' @field storativity aquifer storativity
    storativity = NULL,
    #' @field transmissivity aquifer transmissivity
    transmissivity = NULL,
    #' @field max_terms number of terms to use in Hantush solution.  More is more
    #'  precise but slower.
    max_terms = NULL,
    initialize = function(time,
                          flow_rate,
                          leakage = 100.0,
                          radius = 100.0,
                          storativity = 1e-6,
                          transmissivity = 1e-4,
                          max_terms = 20,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      flow_rate <- deparse(substitute(flow_rate))
      env_list <- get_function_arguments()
      env_list$step_name <- "step_aquifer_leaky"
      env_list$type <- "add"
      super$initialize(
        terms = c(as.symbol(time), as.symbol(flow_rate)),
        env_list
      )


      # step specific values
      self$time <- time
      self$flow_rate <- flow_rate

      self$leakage <- leakage
      self$radius <- radius

      # K & Ss
      self$storativity <- storativity
      self$transmissivity <- transmissivity
      self$max_terms <- max_terms

      self$columns <- c(time, flow_rate)

      invisible(self)
    },
    bake = function(new_data) {
      hantush_jacob(
        new_data[[1]],
        new_data[[2]],
        self$radius,
        self$storativity,
        self$transmissivity,
        self$leakage,
        self$max_terms
      )
    }
  )
)
