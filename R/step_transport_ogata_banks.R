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
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {

      # get function parameters to pass to parent
      step_name    <- "step_transport_ogata_banks"
      type         <- 'add'
      enq <- rlang::enquos(time)
      enq <- rlang::enquos(distance)
      inputs <- c(
        as.list(rlang::quos(...)),
        rlang::env_get_list(env = environment(),
                            formalArgs(super$initialize))
      )
      do.call(super$initialize, inputs)

      # step specific values
      self$time = enquos(time)
      self$distance = enquos(distance)
      self$concentration_initial = concentration_initial
      self$velocity = velocity
      self$diffusion = diffusion
      self$retardation = retardation
      self$decay = decay

      invisible(self)
    },

    bake = function(new_data) {

      ob <- ogata_banks(
        new_data[[1]],
        new_data[[2]],
        self$concentration_initial,
        self$velocity,
        self$diffusion,
        self$retardation,
        self$decay
      )

      names(ob) <- file.path(self$id, new_data[[1]], fsep = '_')

    }
  )
)


