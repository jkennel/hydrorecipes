#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Sudicky and Frind 1982 Parallel Fractures Adapted for Heat Step --------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepTransportFracturesHeat` Sudicky and Frind 1982 solution adapted for
#' heat. Two parallel fractures.
#'
#' @param time vector elapsed time (t)
#' @param distance_fracture vector distance along fracture (z)
#' @param distance_matrix vector distance into matrix (x)
#' @param temperature_influent vector temperature history (t_in)
#' @param time_influent vector time of influent values (t_in)
#' @param temperature_initial double temperature  (t_0)
#' @param fracture_aperture double fracture aperture (2b)
#' @param fracture_spacing double fracture aperture (2B)
#' @param velocity double water velocity in fracture (v)
#' @param thermal_conductivity_water double water thermal conductivity (λ_f)
#' @param thermal_conductivity_solids double solids thermal conductivity (λ_s)
#' @param specific_heat_water double specific heat of water
#' @param specific_heat_solids double specific heat of solid particles
#' @param density_water double density of the water (ρ_w)
#' @param density_solids double density of the solid particles (ρ_s)
#' @param porosity double matrix porosity (θ)
#' @param n_terms integer the number of laplace terms
#'
#' @inheritParams Step
#'
#' @family transport
#'
#' @export
StepTransportFracturesHeat <- R6Class(

  classname = 'step_transport_fractures_heat',
  inherit = Step,

  public = list(

    # step specific variables
    time = NULL,
    distance_fracture = NULL,
    distance_matrix = NULL,
    temperature_influent = NULL,
    time_influent = NULL,
    temperature_initial = NULL,
    fracture_aperture = NULL,
    fracture_spacing = NULL,
    velocity = NULL,
    thermal_conductivity_water = NULL,
    thermal_conductivity_solids = NULL,
    specific_heat_water = NULL,
    specific_heat_solids = NULL,
    density_water = NULL,
    density_solids = NULL,
    porosity = NULL,
    n_terms = NULL,

    initialize = function(time,
                          distance_fracture,
                          distance_matrix,
                          temperature_influent = 15.0,
                          time_influent = 0.0,
                          temperature_initial = 10,
                          fracture_aperture = 2e-4,
                          fracture_spacing = 1.0,
                          velocity = 0.1 / 86400.0,
                          thermal_conductivity_water = 0.615,
                          thermal_conductivity_solids = 3.4,
                          specific_heat_water = 4192,
                          specific_heat_solids = 908,
                          density_water = 1.0,
                          density_solids = 2.5,
                          porosity = 0.1,
                          n_terms = 30L,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      distance_fracture <- deparse(substitute(distance_fracture))
      distance_matrix <- deparse(substitute(distance_matrix))
      env_list <- get_function_arguments()
      env_list$step_name <- 'step_transport_fractures_heat'
      env_list$type <- 'add'
      super$initialize(terms = c(as.symbol(time),
                                 as.symbol(distance_fracture),
                                 as.symbol(distance_matrix)),
                       env_list)

      # step specific values
      self$time = time
      self$distance_fracture = distance_fracture
      self$distance_matrix = distance_matrix
      self$temperature_influent = temperature_influent
      self$time_influent = time_influent
      self$temperature_initial = temperature_initial
      self$fracture_aperture = fracture_aperture
      self$fracture_spacing = fracture_spacing
      self$velocity = velocity
      self$thermal_conductivity_water = thermal_conductivity_water
      self$thermal_conductivity_solids = thermal_conductivity_solids
      self$specific_heat_water = specific_heat_water
      self$specific_heat_solids = specific_heat_solids
      self$density_water = density_water
      self$density_solids = density_solids
      self$porosity = porosity
      self$n_terms = n_terms

      self$columns <- c(time, distance_fracture, distance_matrix)

      invisible(self)
    },

    bake = function(new_data) {

      parallel_fractures_heat(
        new_data[[1]], # time
        new_data[[2]], # z
        new_data[[3]], # x
        self$temperature_influent,
        self$time_influent,
        self$temperature_initial,
        self$fracture_aperture,
        self$fracture_spacing,
        self$velocity,
        self$thermal_conductivity_water,
        self$thermal_conductivity_solids,
        self$specific_heat_water,
        self$specific_heat_solids,
        self$density_water,
        self$density_solids,
        self$porosity,
        self$n_terms)
    }
  )
)


