#' R6 Class
#'
#' `StepTransportOgataBanks` Ogata-Banks solution.
#'
#' @param time vector time
#' @param distance vector x position
#' @param concentration_initial double concentration
#' @param velocity double velocity
#' @param diffusion double diffusion coefficient
#' @param retardation double retardation coefficient
#' @param decay double decay coefficient
#' @inheritParams Step
#'
#' @family transport
#'
#' @export
StepTransportOgataBanks <- R6Class(

  classname = 'step_transport_ogata_banks',
  inherit = Step,

  public = list(

    # step specific variables
    diffusion = NULL,
    retardation = NULL,
    decay = NULL,
    velocity = NULL,
    concentration_initial = NULL,
    distance = NULL,
    time = NULL,

    initialize = function(time,
                          distance,
                          concentration_initial = 1.0,
                          velocity = 0.1,
                          diffusion = 0.1,
                          retardation = 1.0,
                          decay = 0.0,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      distance <- deparse(substitute(distance))
      env_list <- get_function_arguments()
      env_list$step_name <- 'step_transport_ogata_banks'
      env_list$type <- 'add'
      super$initialize(terms = c(as.symbol(time), as.symbol(distance)),
                       env_list)

      # step specific values
      self$time = time
      self$distance = distance
      self$concentration_initial = concentration_initial
      self$velocity = velocity
      self$diffusion = diffusion
      self$retardation = retardation
      self$decay = decay

      self$columns <- c(time, distance)

      invisible(self)
    },

    bake = function(new_data) {

      ogata_banks_list(
        new_data[[1]],
        new_data[[2]],
        self$concentration_initial,
        self$velocity,
        self$diffusion,
        self$retardation,
        self$decay
      )

    }
  )
)


