#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Generate Synthetic Earth Tides Step ------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' R6 Class
#'
#' `StepEarthtide` generates sin and cosine wavegroups.
#'
#' @inheritParams Step, earthtide::calc_earthtide
#'
#' @export
StepEarthtide <- R6Class(
  classname = "step_earthtide",
  inherit = Step,
  public = list(

    # step specific variables
    latitude = NA_real_,
    longitude = NA_real_,
    elevation = NA_real_,
    azimuth = NA_real_,
    gravity = NA_real_,
    earth_radius = NA_real_,
    earth_eccentricity = NA_real_,
    cutoff = NA_real_,
    catalog = NA_character_,
    eop = NULL,
    scale = NA,
    n_thread = NA_integer_,
    do_predict = NA,
    method = NA_character_,
    initialize = function(terms,
                          do_predict = TRUE,
                          method = "gravity",
                          latitude = 0.0,
                          longitude = 0.0,
                          elevation = 0.0,
                          azimuth = 0.0,
                          gravity = 0.0,
                          earth_radius = 6378136.3,
                          earth_eccentricity = 0.0066943979514,
                          cutoff = 1e-6,
                          catalog = "ksm04",
                          eop = NULL,
                          scale = TRUE,
                          n_thread = 1L,
                          role = "predictor",
                          ...) {
      # get function parameters to pass to parent
      terms <- substitute(terms)
      env_list <- get_function_arguments()
      env_list$step_name <- "step_earthtide"
      env_list$type <- "add"
      super$initialize(
        terms = terms,
        env_list[names(env_list) != "terms"]
      )

      # step specific values
      self$method <- method
      self$latitude <- latitude
      self$longitude <- longitude
      self$elevation <- elevation
      self$azimuth <- azimuth
      self$gravity <- gravity
      self$earth_radius <- earth_radius
      self$earth_eccentricity <- earth_eccentricity
      self$cutoff <- cutoff
      self$catalog <- catalog
      self$eop <- eop
      self$scale <- scale
      self$n_thread <- n_thread
      self$do_predict <- do_predict


      invisible(self)
    },
    bake = function(new_data) {
      column_name <- self$columns

      et <- calc_earthtide(
        new_data[[column_name]],
        do_predict = self$do_predict,
        method = self$method,
        latitude = self$latitude,
        longitude = self$longitude,
        elevation = self$elevation,
        azimuth = self$azimuth,
        gravity = self$gravity,
        earth_radius = self$earth_radius,
        earth_eccentricity = self$earth_eccentricity,
        cutoff = self$cutoff,
        catalog = self$catalog,
        eop = self$eop,
        scale = self$scale,
        n_thread = self$n_thread
      )

      names(et) <- file.path(self$id,
        names(et),
        fsep = "_"
      )

      et
    }
  )
)
