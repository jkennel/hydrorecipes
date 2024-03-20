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
    wave_groups = NA_real_,
    catalog = NA_character_,
    eop = NULL,
    scale = NA,
    n_thread = NA_integer_,
    do_predict = NA,
    method = NA_character_,
    frequency = NA_real_,
    return_matrix = TRUE,

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
                          wave_groups = data.frame(start = 0, end = 8),
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
        env_list[names(env_list) != "terms"],
        ...
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
      self$wave_groups <- wave_groups
      self$catalog <- catalog
      self$eop <- eop
      self$scale <- scale
      self$n_thread <- n_thread
      self$do_predict <- do_predict

      if (!do_predict) {

        self$frequency = earthtide::get_main_frequency(wave_groups[[1]],
                                                       wave_groups[[2]])

      } else {
        self$frequency = NA_real_
      }


      invisible(self)
    },
    bake = function(new_data) {
      column_name <- self$columns
      et <- mctl(calc_earthtide(
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
        wave_groups = self$wave_groups,
        catalog = self$catalog,
        eop = self$eop,
        scale = self$scale,
        return_matrix = self$return_matrix,
        n_thread = self$n_thread
      ))

      if (self$do_predict) {
        self$new_columns <- paste0(self$prefix)
      } else {
        self$new_columns <- paste(rep(self$prefix, length(self$frequency) * 2L),
               rep(c("cos", "sin"), length(self$frequency)),
               rep(1:length(self$frequency), each = 2L), sep = "_")
      }

      names(et) <- self$new_columns

      et
    },
    response = function(co) {

      if (is.na(self$frequency)) {
        super$response()
      }

      f <- self$frequency
      n <- length(f)
      x <- rep(f, 2)

      cos_coefficient <- seq(1, n, 2)
      sin_coefficient <- seq(2, n, 2)
      amp_phase <- c(
        sqrt(cos_coefficient^2 + sin_coefficient^2), # amplitude
        atan2(cos_coefficient, sin_coefficient)      # phase
      )
      variable <- c(
        rep("amplitude", n),
        rep("phase", n)
      )

      list(x, variable, value = amp_phase, step_id = self$id)
    }
  )
)
