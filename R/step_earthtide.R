#' Earth tide response
#'
#' `step_earthtide` creates a *specification* of a recipe step
#'  that are the Earth tide harmonics for a particular
#'  location. This step requires the [earthtide package](https://CRAN.R-project.org/package=earthtide).
#'
#' @inheritParams recipes::step_lag
#' @inheritParams earthtide::calc_earthtide
#' @param ... One or more selector functions to choose which
#'  variables are affected by the step. See [selections()]
#'  for more details. For the `tidy` method, these are not
#'  currently used.
#' @param role Defaults to "earthtide"
#' @param prefix A prefix for generated column names, default to
#'  "earthtide_".
#' @param keep_original_cols A logical to keep the original variables in the
#'  output. Defaults to `FALSE`.
#' @return An updated version of `recipe` with the new step
#'  added to the sequence of existing steps (if any).
#' @keywords datagen
#' @concept generate Earth tide harmonics
#' @export
#'
#'
#' @details There are many waves (thousands) that make up a tidal signal.
#'   `step_earthtide` calculates the Earth tide signal for a time and
#'    location. The tidal signal can be estimated as a single summed curve when
#'    `do_predict = TRUE` or as a set of wave groups when `do_predict = FALSE`.
#'    Wave groups are ranges of frequencies identified by start and end
#'    frequencies.  For example, if you have one month of data the M2 signal
#'    could be described as the sum of all the waves between 1.914129 and
#'    1.950419 cycles per day. The regressors for each wave group have a sin
#'    and cos component. See [recipes::step_harmonic()] for a simplified
#'    version of this where each sin and cos curve corresponds to a single wave.
#'
#' @examples
#' library(earthtide)
#' data(eterna_wavegroups)
#' data(transducer)
#'
#' transducer <- transducer[, c('datetime', 'wl'),]
#' t_sub <- transducer[(as.numeric(transducer$datetime) %% 14400) == 0, ]
#' wg <- na.omit(eterna_wavegroups[eterna_wavegroups$time == '1 month',])
#'
#' recipe(wl ~ ., data = t_sub) |>
#'   step_earthtide(datetime,
#'                  latitude = 34,
#'                  longitude = -118.5,
#'                  wave_groups = wg,
#'                  do_predict = FALSE) |>
#'   prep()
#'
#' recipe(wl ~ ., data = t_sub) |>
#'   step_earthtide(datetime,
#'                  latitude = 34,
#'                  longitude = -118.5,
#'                  wave_groups = wg,
#'                  do_predict = TRUE) |>
#'   prep()
#'
#' @seealso [step_earthtide()] [earthtide::calc_earthtide()]
#'  [recipes::step_harmonic()]
#' @importFrom recipes add_step step recipes_eval_select ellipse_check rand_id
step_earthtide <-
  function(recipe,
           ...,
           role = "predictor",
           trained = FALSE,
           do_predict = FALSE,
           method = "gravity",
           astro_update = 1L,
           latitude = 0,
           longitude = 0,
           elevation = 0,
           azimuth = 0,
           gravity = 0,
           earth_radius = 6378136.3,
           earth_eccen = 0.0066943979514,
           cutoff = 1e-06,
           wave_groups = NULL,
           catalog = "ksm04",
           eop = NULL,
           scale = TRUE,
           prefix = "earthtide_",
           columns = NULL,
           keep_original_cols = FALSE,
           skip = FALSE,
           id = rand_id("earthtide")) {
    add_step(
      recipe,
      step_earthtide_new(
        terms = ellipse_check(...),
        role = role,
        trained = trained,
        do_predict = do_predict,
        method = method,
        astro_update = astro_update,
        latitude = latitude,
        longitude = longitude,
        elevation = elevation,
        azimuth = azimuth,
        gravity = gravity,
        earth_radius = earth_radius,
        earth_eccen = earth_eccen,
        cutoff = cutoff,
        wave_groups = wave_groups,
        catalog = catalog,
        eop = eop,
        scale = scale,
        prefix = prefix,
        columns = columns,
        keep_original_cols = keep_original_cols,
        skip = skip,
        id = id
      )
    )
  }

step_earthtide_new <-
  function(terms, role, trained, do_predict, method, astro_update, latitude, longitude,
           elevation, azimuth, gravity, earth_radius, earth_eccen, cutoff,
           wave_groups, catalog, eop, scale, prefix, columns,
           keep_original_cols,skip, id) {
    step(
      subclass = "earthtide",
      terms = terms,
      role = role,
      trained = trained,
      do_predict = do_predict,
      method = method,
      astro_update = astro_update,
      latitude = latitude,
      longitude = longitude,
      elevation = elevation,
      azimuth = azimuth,
      gravity = gravity,
      earth_radius = earth_radius,
      earth_eccen = earth_eccen,
      cutoff = cutoff,
      wave_groups = wave_groups,
      catalog = catalog,
      eop = eop,
      scale = scale,
      prefix = prefix,
      columns = columns,
      keep_original_cols = keep_original_cols,
      skip = skip,
      id = id
    )
  }

#' @export
prep.step_earthtide <- function(x, training, info = NULL, ...) {

  step_earthtide_new(
    terms = x$terms,
    role = x$role,
    trained = TRUE,
    do_predict = x$do_predict,
    method = x$method,
    astro_update = x$astro_update,
    latitude = x$latitude,
    longitude = x$longitude,
    elevation = x$elevation,
    azimuth = x$azimuth,
    gravity = x$gravity,
    earth_radius = x$earth_radius,
    earth_eccen = x$earth_eccen,
    cutoff = x$cutoff,
    wave_groups = x$wave_groups,
    catalog = x$catalog,
    eop = x$eop,
    scale = x$scale,
    prefix = x$prefix,
    columns = recipes_eval_select(x$terms, training, info = info),
    keep_original_cols = get_keep_original_cols(x),
    skip = x$skip,
    id = x$id
  )

}

#' @importFrom dplyr bind_cols
#' @importFrom tibble as_tibble
#' @importFrom earthtide calc_earthtide
#' @importFrom recipes bake prep
#'
#' @export
bake.step_earthtide <- function(object, new_data, ...) {

  et <- as_tibble(calc_earthtide(utc = new_data[[object$columns]],
                                 do_predict = object$do_predict,
                                 method = object$method,
                                 astro_update = object$astro_update,
                                 latitude = object$latitude,
                                 longitude = object$longitude,
                                 elevation = object$elevation,
                                 azimuth = object$azimuth,
                                 gravity = object$gravity,
                                 earth_radius = object$earth_radius,
                                 earth_eccen = object$earth_eccen,
                                 cutoff = object$cutoff,
                                 wave_groups = object$wave_groups,
                                 catalog = object$catalog,
                                 eop = object$eop,
                                 scale = object$scale,
                                 return_matrix = TRUE))
  if(object$do_predict) {
    names(et) <- paste0(object$prefix, object$columns)
  } else {
    names(et) <- paste0(object$prefix,
                        rep(c('cos', 'sin'), times = nrow(object$wave_groups)),
                        '_', rep(1:nrow(object$wave_groups), each = 2))
  }

  bind_cols(new_data, et)

}

#' @importFrom recipes printer
print.step_earthtide <-
  function(x, width = max(20, options()$width - 30), ...) {
    cat("earthtide ",  sep = "")
    printer(x$columns, x$terms, x$trained, width = width)
    invisible(x)
  }


#' @export
tidy.step_earthtide <- function(x, ...) {
 tidy2(x, ...)
}

#' @rdname tidy2.recipe
#' @export
tidy2.step_earthtide <- function(x, ...) {

  n_terms <- nrow(x$wave_groups)
  res <-
    tibble(terms = rep(sel2char(x$terms), n_terms * 2),
           latitude = rep(x$latitude, n_terms * 2),
           longitude = rep(x$longitude, n_terms * 2),
           elevation = rep(x$elevation, n_terms * 2),
           azimuth = rep(x$azimuth, n_terms * 2),
           gravity = rep(x$gravity, n_terms * 2),
           cutoff = rep(x$cutoff, n_terms * 2),
           catalog = rep(x$catalog, n_terms * 2),
           frequency_start = rep(x$wave_groups$start, times = 2),
           frequency_end = rep(x$wave_groups$end, times = 2),
           frequency = rep(earthtide::get_main_frequency(x$wave_groups$start,
                                                         x$wave_groups$end), times = 2)
    )

  res$key <- paste0(x$prefix,
                    rep(c('sin_', 'cos_'), each = n_terms),
                    rep(1:n_terms, times = 2))

  res$id <- x$id
  res$step_name <- 'step_earthtide'

  res
}

#' @export
tidy.step_earthtide <- function(x, ...) {

  tidy2.step_earthtide(x, ...)
}
