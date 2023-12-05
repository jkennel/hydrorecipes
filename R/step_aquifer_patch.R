#' R6 Class
#'
#' `StepAquiferPatch` barker_herbert.
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
#' @inheritParams Step
#' @family aquifer
#' @export
StepAquiferPatch <- R6Class(

  classname = 'step_aquifer_patch',
  inherit = Step,

  public = list(

    # step specific variables
    time = NULL,
    flow_rate = NULL,
    thickness = NULL,
    radius = NULL,
    radius_patch = NULL,
    specific_storage_inner = NULL,
    specific_storage_outer = NULL,
    hydraulic_conductivity_inner = NULL,
    hydraulic_conductivity_outer = NULL,

    # calculated
    storativity_inner = NULL,
    storativity_outer = NULL,
    transmissivity_inner = NULL,
    transmissivity_outer = NULL,
    n_stehfest = NULL,

    initialize = function(time,
                          flow_rate = 0.01,
                          thickness = 1.0,
                          radius = 200.0,
                          radius_patch = 100.0,
                          specific_storage_inner = 1.0e-6,
                          specific_storage_outer = 1.0e-5,
                          hydraulic_conductivity_inner = 1.0e-4,
                          hydraulic_conductivity_outer = 1.0e-6,
                          n_stehfest = 12L,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {

      # get function parameters to pass to parent
      step_name    <- "step_aquifer_patch"
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
      self$flow_rate = flow_rate

      self$thickness = thickness
      self$radius = radius
      self$radius_patch = radius_patch

      # K & Ss
      self$specific_storage_inner = specific_storage_inner
      self$specific_storage_outer = specific_storage_outer
      self$hydraulic_conductivity_inner = hydraulic_conductivity_inner
      self$hydraulic_conductivity_outer = hydraulic_conductivity_outer

      # T & S
      self$storativity_inner = specific_storage_inner * thickness
      self$storativity_outer = specific_storage_outer * thickness
      self$transmissivity_inner = hydraulic_conductivity_inner * thickness
      self$transmissivity_outer = hydraulic_conductivity_outer * thickness

      self$n_stehfest = n_stehfest


      invisible(self)
    },
#(frecipes:::stehfest_barker_herbert(as.numeric(1:n),1.0, 100.0, 200.0, 1e-3, 1e-3, 1e-5, 1e-5, 12L)[[1]])

    bake = function(new_data) {
      setNames(stehfest_barker_herbert(
          new_data[[1]],
          self$flow_rate,
          self$radius,
          self$radius_patch,
          self$transmissivity_inner,
          self$transmissivity_outer,
          self$storativity_inner,
          self$storativity_outer,
          self$n_stehfest
        ),
        self$id)
    }
  )
)




