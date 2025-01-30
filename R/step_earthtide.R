#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# Generate Synthetic Earth Tides Step ------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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
    earth_eccen = NA_real_,
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
    astro_update = NA_integer_,
    interp_factor = NA_integer_,
    utc_interp = NULL,

    initialize = function(terms,
                          do_predict = TRUE,
                          method = "gravity",
                          latitude = 0.0,
                          longitude = 0.0,
                          elevation = 0.0,
                          azimuth = 0.0,
                          gravity = 0.0,
                          earth_radius = 6378136.3,
                          earth_eccen = 0.0066943979514,
                          cutoff = 1e-6,
                          wave_groups = data.frame(start = 0, end = 8),
                          catalog = "ksm04",
                          eop = NULL,
                          scale = TRUE,
                          n_thread = 1L,
                          astro_update = 1L,
                          interp_factor = 1L,
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
      self$earth_eccen <- earth_eccen
      self$cutoff <- cutoff
      self$wave_groups <- wave_groups
      self$catalog <- catalog
      self$eop <- eop
      self$scale <- scale
      self$n_thread <- n_thread
      self$do_predict <- do_predict
      self$astro_update <- astro_update
      self$interp_factor = interp_factor

      if (!do_predict) {

        self$frequency = earthtide::get_main_frequency(wave_groups[[1L]],
                                                       wave_groups[[2L]])

      } else {
        self$frequency = NA_real_
      }


      invisible(self)
    },
    bake = function(s) {

      column_name <- self$columns

      if (self$interp_factor != 1L) {
        self$utc_interp <- s[["result"]][[self$columns]]
        utc <- self$utc_interp[(0:(length(self$utc_interp) - 1L) %% self$interp_factor) == 0L]
      } else {
        utc <- s[["result"]][[self$columns]]
      }



      et <- unclass(earthtide::calc_earthtide(
        utc = utc,
        do_predict = self$do_predict,
        method = self$method,
        latitude = self$latitude,
        longitude = self$longitude,
        elevation = self$elevation,
        azimuth = self$azimuth,
        gravity = self$gravity,
        earth_radius = self$earth_radius,
        earth_eccen = self$earth_eccen,
        cutoff = self$cutoff,
        wave_groups = self$wave_groups,
        catalog = self$catalog,
        eop = self$eop,
        scale = self$scale,
        return_matrix = FALSE,
        n_thread = self$n_thread,
        astro_update = self$astro_update,
        utc_interp = self$utc_interp
      )[, -1L, drop = FALSE])



      if (self$do_predict) {
        self$new_columns <- paste0(self$prefix)
      } else {
        self$new_columns <- paste(rep(self$prefix, length(self$frequency) * 2L),
               rep(c("cos", "sin"), length(self$frequency)),
               rep(1:length(self$frequency), each = 2L), sep = "_")
      }

      names(et) <- self$new_columns

      self$result <- et
      return(NULL)

    },
    response = function(co) {

      if (all(is.na(self$frequency))) {
        return(super$response(co))
      }

      f  <- self$frequency

      nr <- length(f)
      nc <- ncol(co)

      x  <- rep(rep(f, each = nc), 2L)

      cos_coefficient <- co[seq.int(1L, nr * 2L, 2L), , drop = FALSE]
      sin_coefficient <- co[seq.int(2L, nr * 2L, 2L), , drop = FALSE]

      amp_phase <- c(
        as.vector(t(sqrt(cos_coefficient^2 + sin_coefficient^2))), # amplitude
        as.vector(t(atan2(cos_coefficient, sin_coefficient)))      # phase
      )

      variable <- c(
        rep.int("amplitude", nr * nc),
        rep.int("phase", nr * nc)
      )

      list(x = x,
           variable = variable,
           value = amp_phase,
           step_id = rep.int(self$id, nr * nc),
           outcome = rep(colnames(co), each = nr),
           term = rep.int(self$prefix, nr * nc),
           step_columns = rep.int(paste(self$columns, collapse = "_"), nr * nc))

    }
  )
)
