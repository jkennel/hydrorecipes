#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Sudicky and Frind 1982 Parallel Fractures Step -------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
StepTransportFracturesSolute <- R6Class(

  classname = 'step_transport_fractures_solute',
  inherit = Step,

  public = list(

    # step specific variables
    time = NULL,
    distance_fracture = NULL,
    distance_matrix = NULL,
    concentration_influent = NULL,
    time_influent = NULL,
    concentration_initial = NULL,
    fracture_aperture = NULL,
    fracture_spacing = NULL,
    velocity = NULL,
    dispersivity_longitudinal = NULL,
    diffusion = NULL,
    sorption_fracture = NULL,
    sorption_matrix = NULL,
    decay = NULL,
    density_bulk = NULL,
    porosity = NULL,
    tortuosity = NULL,
    n_terms = NULL,

    initialize = function(time,
                          distance_fracture,
                          distance_matrix,
                          concentration_influent = 1.0,
                          time_influent = 0.0,
                          concentration_initial = 0.0,
                          fracture_aperture = 2e-4,
                          fracture_spacing = 1.0,
                          velocity = 0.1 / 86400.0,
                          dispersivity_longitudinal = 0.1,
                          diffusion = 1e-9,
                          sorption_fracture = 0.0,
                          sorption_matrix = 0.0,
                          decay = 1e15, # no decay
                          density_bulk = 2.5,
                          porosity = 0.10,
                          tortuosity = 0.1,
                          n_terms = 30L,
                          role = "predictor",
                          ...) {

      # get function parameters to pass to parent
      time <- deparse(substitute(time))
      distance_fracture <- deparse(substitute(distance_fracture))
      distance_matrix <- deparse(substitute(distance_matrix))
      env_list <- get_function_arguments()
      env_list$step_name <- 'step_transport_fractures_solute'
      env_list$type <- 'add'
      super$initialize(terms = c(as.symbol(time),
                                 as.symbol(distance_fracture),
                                 as.symbol(distance_matrix)),
                       env_list,
                       ...)

      # step specific values
      self$time = time
      self$distance_fracture = distance_fracture
      self$distance_matrix = distance_matrix
      self$concentration_influent = concentration_influent
      self$time_influent = time_influent
      self$concentration_initial = concentration_initial
      self$fracture_aperture = fracture_aperture
      self$fracture_spacing = fracture_spacing
      self$velocity = velocity
      self$dispersivity_longitudinal = dispersivity_longitudinal
      self$diffusion = diffusion
      self$sorption_fracture = sorption_fracture
      self$sorption_matrix = sorption_matrix
      self$decay = decay
      self$density_bulk = density_bulk
      self$porosity = porosity
      self$tortuosity = tortuosity
      self$n_terms = n_terms

      self$columns <- c(time, distance_fracture, distance_matrix)

      invisible(self)
    },

    bake = function(new_data) {

      pfs <- list(parallel_fractures_solute(
        new_data[[1]], # time
        new_data[[2]], # z
        new_data[[3]], # x
        self$concentration_influent,
        self$time_influent,
        self$concentration_initial,
        self$fracture_aperture,
        self$fracture_spacing,
        self$velocity,
        self$dispersivity_longitudinal,
        self$diffusion,
        self$sorption_fracture,
        self$sorption_matrix,
        self$decay,
        self$density_bulk,
        self$porosity,
        self$tortuosity,
        self$n_terms))


      self$new_columns <- self$prefix

      names(pfs) <- self$new_columns

      self$result <- pfs
      self$result

    }
  )
)


