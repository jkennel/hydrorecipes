#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# recipe -----------------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Create a new R6 recipe. This is analogous to the the list structure that the
#' *recipes* package uses.
#'
#' @inheritParams stats::lm
#' @param ... additional arguments to pass to Recipe$new().  This is currently
#' not used.
#'
#' @return R6 recipe object
#' @export
#'
#' @examples
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat)
#'
recipe <- function(formula, data, ...) {
  Recipe$new(formula, data, ...)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# steps ------------------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_add_vars
#'
#' @description
#'   Add a variable from the initial data set after recipe creation.
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#' dat <- data.frame(x = rnorm(10), y = rnorm(10), z = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_add_vars(z) |> plate()
#'
#' rec <- recipe(y~x, data = dat) |>
#'        plate()
step_add_vars <- function(.rec,
                          terms,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAddVars$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_add_noise
#'
#' @description
#'   Add noise.
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
step_add_noise <- function(.rec,
                          terms,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAddNoise$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_aquifer_constant_drawdown
#'
#' @description
#'   Jacob and Lohman solution for constant drawdown test.
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#' dat <- data.frame(x = rnorm(10), y = rnorm(10), z = rnorm(10))
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_aquifer_constant_drawdown(time = times,
#'                                  drawdown = 10,
#'                                  thickness = 10,
#'                                  radius_well = 0.15,
#'                                  specific_storage = 1e-6,
#'                                  hydraulic_conductivity = 1,
#'                                  n_terms = 12L)
#'
step_aquifer_constant_drawdown <- function(.rec,
                                           time,
                                           drawdown = 1.0,
                                           thickness = 1.0,
                                           radius_well = 0.15,
                                           specific_storage = 1.0e-6,
                                           hydraulic_conductivity = 1.0e-4,
                                           n_terms = 16,
                                           role = "predictor",
                                           ...) {
  time <- substitute(time)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAquiferConstantDrawdown$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_aquifer_grf
#'
#' @description
#' Generates the drawdown using the Generalized
#' Radial Flow (GRF) model. This method defaults to a fast FFT
#' convolution so many rates can be included, but requires a regular
#' time series.
#'
#' @inheritParams step_scale
#'
#' @return The drawdown using the GRF model
#'
#' @references
#' Barker, J.A., A generalized radial flow model for hydraulic tests
#'  in fractured rock. Water Resour. Res., 24 (1988), pp. 1796-1804,
#'  10.1029/WR024i010p01796
#'
#' @family aquifer
#'
#' @examples
#' dat <- data.frame(x = as.numeric(1:rows),
#'                   y = rep(0.01, rows))
#' formula <- as.formula(y~x)
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_aquifer_grf(time = x, flow_rate = y)
#'
#' @export
step_aquifer_grf <- function(.rec,
                             time,
                             flow_rate,
                             thickness = 1.0,
                             radius = 100.0,
                             specific_storage = 1.0e-6,
                             hydraulic_conductivity = 1.0e-4,
                             flow_dimension = 2.0,
                             role = "predictor",
                             ...) {
  time <- substitute(time)
  flow_rate <- substitute(flow_rate)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAquiferGRF$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_aquifer_theis
#'
#' @description
#; Generates the drawdown using the Generalized
#' Radial Flow (GRF) model with flow dimension equal to 2. This method defaults
#' to a fast FFT convolution so many rates can be included, but requires a
#' regular time series.
#'
#' @inheritParams step_scale
#' @inheritParams step_aquifer_grf
#'
#' @return The drawdown using the Theis model
#'
#' @references
#' Barker, J.A., A generalized radial flow model for hydraulic tests
#'  in fractured rock. Water Resour. Res., 24 (1988), pp. 1796-1804,
#'  10.1029/WR024i010p01796
#'
#' Theis, C.V., 1935: The relation between the lowering of the piezometric
#'  surface and the rate and duration of discharge of a well using
#'  ground-water storage, Transactions of the American Geophysical Union,
#'  16th Annual Meeting, Part 2, pp. 519-524.
#'
#' @family aquifer
#'
#' @examples
#' dat <- data.frame(x = as.numeric(1:rows),
#'                   y = rep(0.01, rows))
#' formula <- as.formula(y~x)
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_aquifer_theis(time = x, flow_rate = y)
#'
#' @export
step_aquifer_theis <- function(.rec,
                               time,
                               flow_rate,
                               thickness = 1.0,
                               radius = 100.0,
                               specific_storage = 1.0e-6,
                               hydraulic_conductivity = 1.0e-4,
                               flow_dimension = 2.0,
                               role = "predictor",
                               ...) {
  time <- substitute(time)
  flow_rate <- substitute(flow_rate)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAquiferTheis$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_aquifer_leaky
#'
#' @description
#' hantush_jacob 1955 solution for a leaky aquifer
#'
#' @inheritParams step_scale
#' @inheritParams step_aquifer_grf
#'
#' @return The drawdown using the Hantush and Jacob 1955 model
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
#' @examples
#' dat <- data.frame(x = as.numeric(1:rows),
#'                   y = rep(0.01, rows))
#' formula <- as.formula(y~x)
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_aquifer_leaky(time = x, flow_rate = y)
#'
#' @export
step_aquifer_leaky <- function(.rec,
                               time,
                               flow_rate,
                               leakage = 100.0,
                               radius = 100.0,
                               storativity = 1e-6,
                               transmissivity = 1e-4,
                               max_terms = 20,
                               role = "predictor",
                               ...) {
  time <- substitute(time)
  flow_rate <- substitute(flow_rate)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAquiferLeaky$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_baro_acworth
#'
#' @description
#' Acworth 2016 frequency based method for calculating barometric efficiency
#' in the presence of Earth tides
#'
#' @inheritParams step_fft_pgram
#'
#' @param water_level \code{variable} unquoted water level column name
#' @param barometric_pressure \code{variable} unquoted barometric pressure
#'   column name
#' @param earth_tide \code{variable} unquoted Earth tide column name
#' @param frequency_a \code{double} Earth tide frequency
#' @param frequency_b \code{double} Related barometric frequency
#' @param inverse \code{logical} whether the barometric relationship is inverse
#'
#' @return \code{double} barometric efficiency using Acworth's method
#'
#' @family barometric
#'
#' @references
#' Acworth, R.I., Halloran, L.J., Rau, G.C., Cuthbert, M.O. and Bernardi, T.L.,
#'  2016. An objective frequency domain method for quantifying confined aquifer
#'  compressible storage using Earth and atmospheric tides. Geophysical Research
#'  Letters, 43(22), pp.11-671.
#'
#' @examples
#'
#' @export
step_baro_acworth <- function(.rec,
                              water_level,
                              barometric_pressure,
                              earth_tide,
                              frequency_a = 1.9324, # m2
                              frequency_b = 2.0,    # s2
                              inverse = FALSE,
                              spans = 5,
                              detrend = TRUE,
                              demean = TRUE,
                              taper = 0.1,
                              role = "augment",
                              ...) {
  water_level <- substitute(water_level)
  barometric_pressure <- substitute(barometric_pressure)
  earth_tide <- substitute(earth_tide)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepBaroAcworth$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_baro_clark
#'
#' @description
#' Clark 1967 solution for calculating barometric efficiency (Algorithm from Batu 1998, pg 76)
#'
#' @inheritParams step_scale
#'
#' @param dep \code{numeric vector} of the dependent variable (ie:water level)
#' @param ind \code{numeric vector} of the independent variable (ie:barometric pressure)
#' @param lag_space \code{integer} spacing for lags, useful for higher frequency monitoring
#' @param inverse \code{logical} whether the barometric relationship is inverse
#'
#' @return barometric efficiency using Clark's method
#'
#' @family barometric
#'
#' @references
#' Clark, W.E., 1967. Computing the barometric efficiency of a well. Journal
#' of the Hydraulics Division, 93(4), pp.93-98.
#'
#' Batu, V., 1998. Aquifer hydraulics: a comprehensive guide to hydrogeologic
#' data analysis. John Wiley & Sons.
#'
#' @examples
#' dat <- data.frame(x = as.numeric(1:rows),
#'                   y = rep(0.01, rows))
#' formula <- as.formula(y~x)
#'
#' @export
step_baro_clark <- function(.rec,
                            water_level,
                            barometric_pressure,
                            lag_space = 1L,
                            inverse = FALSE,
                            role = "augment",
                            ...) {
  water_level <- substitute(water_level)
  barometric_pressure <- substitute(barometric_pressure)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepBaroClark$new,
                        env_list))
}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_baro_harmonic
#'
#' @description
#'
#' @inheritParams step_fft_pgram
#' @inheritParams step_harmonic
#'
#' @param water_level \code{variable} unquoted water level column name
#' @param barometric_pressure \code{variable} unquoted barometric pressure
#'   column name
#' @param earth_tide \code{variable} unquoted Earth tide column name
#' @param inverse \code{logical} whether the barometric relationship is inverse
#'
#' @return \code{double} barometric efficiency using Acworth's method
#'
#' @family barometric
#'
#' @references
#' Acworth, R.I., Halloran, L.J., Rau, G.C., Cuthbert, M.O. and Bernardi, T.L.,
#'  2016. An objective frequency domain method for quantifying confined aquifer
#'  compressible storage using Earth and atmospheric tides. Geophysical Research
#'  Letters, 43(22), pp.11-671.
#'
#' @examples
#'
#' @export
step_baro_harmonic <- function(.rec,
                               time,
                               water_level,
                               barometric_pressure,
                               earth_tide,
                               frequency = c(1.9324, 2.0),
                               cycle_size = 86400,
                               start = 0.0,
                               inverse = TRUE,
                               role = "augment",
                               ...) {
  time <- substitute(time)
  water_level <- substitute(water_level)
  barometric_pressure <- substitute(barometric_pressure)
  earth_tide <- substitute(earth_tide)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepBaroHarmonic$new,
                        env_list))
}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_aquifer_patch
#'
#' @description
#' barker_herbert 1982 solution for radial patches.
#'
#' @inheritParams step_scale
#' @inheritParams step_aquifer_grf
#'
#' @return The drawdown using the Theis model
#'
#' @references
#' Barker, J.A., and R. Herbert, 1982: Pumping tests in patchy
#'  aquifers, Ground Water, vol. 20, No. 2, pp. 150-155.
#'
#' Butler, J.J., 1988: Pumping tests in nonuniform aquifers – The radially
#'  symmetric case, Journal of Hydrology, Vol. 101, pp. 15-30.
#'
#' @family aquifer
#'
#' @examples
#' dat <- data.frame(x = as.numeric(1:rows),
#'                   y = rep(0.01, rows))
#' formula <- as.formula(y~x)
#'
#' @export
step_aquifer_patch <- function(.rec,
                               time,
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
                               ...) {
  time <- substitute(time)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAquiferPatch$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_center
#' @description
#'   Adds a step to center a data column(s)
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_center(x)
#'
step_center <- function(.rec,
                        terms,
                        role = "predictor",
                        skip = FALSE,
                        na_rm = TRUE,
                        fun = collapse::fmean,
                        keep_original_cols = FALSE,
                        ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepCenter$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_check_na
#'
#' @description
#'   Check columns for NA
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#'
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_scale(x)
#'
step_check_na <- function(.rec,
                          terms,
                          role = "check",
                          ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepCheckNA$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_check_spacing
#'
#' @description
#'   Check columns for NA
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#'
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_scale(x)
#'
step_check_spacing <- function(.rec,
                               terms,
                               role = "check",
                               ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepCheckSpacing$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_distributed_lag
#'
#' @description
#'   generates distributed lag vectors.
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
step_distributed_lag <- function(.rec,
                                 terms,
                                 n_lag = 12L,
                                 max_lag = 86400L,
                                 knots = NA_real_,
                                 basis_matrix = NA_real_,
                                 intercept = FALSE,
                                 role = "predictor",
                                 ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepDistributedLag$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_drop_columns
#'
#' @description
#'   generates distributed lag vectors.
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
step_drop_columns <- function(.rec,
                              terms,
                              role = "modify",
                              ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepDropColumns$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_dummy
#'
#' @description
#'   dummy encode a factor or factor like variable.
#'
#'
#' @inheritParams step_scale
#' @param one_hot logical - use one hot encoding.
#'
#' @return
#' @export
#'
#' @examples
#' dat <- data.frame(x = qF(sample(1:10, 100, replace = TRUE)),
#'                   y = rnorm(100))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_dummy(x, one_hot = FALSE)
#' rec <- recipe(y~x, data = dat) |>
#'        step_dummy(x, one_hot = TRUE)
step_dummy <- function(.rec,
                       terms,
                       one_hot = FALSE,
                       role = "predictor",
                       skip = FALSE,
                       keep_original_cols = FALSE,
                       ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepDummy$new,
                        env_list))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_earthtide
#'
#' @description
#'   Generate synthetic Earth tide waves and wave groups.
#'
#' @inheritParams step_scale, earthtide::calc_earthtide
#'
#' @return
#' @export
#'
#' @examples
step_earthtide <- function(.rec,
                           terms,
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
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepEarthtide$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_fft_coherence
#'
#' @description
#'   estimates the coherence between terms.
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
step_fft_coherence <- function(.rec,
                               terms,
                               role = "predictor",
                               ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepCoherence$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_fft_pgram
#'
#' @description
#'   Transfer function using pgram method.
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
step_fft_pgram <- function(.rec,
                           terms,
                           spans = 3,
                           detrend = TRUE,
                           demean = TRUE,
                           lst = TRUE,
                           taper = 0.1,
                           role = "predictor",
                           ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepPgram$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_fft_welch
#'
#' @description
#'  calculates the periodogram (estimate of spectral density) using
#'  Welch's method.
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
step_fft_welch <- function(.rec,
                           terms,
                           length_subset,
                           overlap = 0.8,
                           window,
                           role = "predictor",
                           ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepWelch$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_fft_transfer_pgram
#'
#' @description
#'  calculates the transfer function using pgram method.
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
step_fft_transfer_pgram <- function(.rec,
                                    terms,
                                    spans = 3,
                                    detrend = TRUE,
                                    demean = TRUE,
                                    taper = 0.1,
                                    role = "predictor",
                                    ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepTransferPgram$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_fft_transfer_welch
#'
#' @description
#'  calculates the transfer function using Welch's method.
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
step_fft_transfer_welch <- function(.rec,
                                    terms,
                                    length_subset,
                                    overlap = 0.8,
                                    window,
                                    role = "predictor",
                                    ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepTransferWelch$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_find_interval
#'
#' @description
#' divides a series into intervals and then performs dummy encoding.
#'
#' @param vec a vector of break points#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
step_find_interval <- function(.rec,
                               terms,
                               vec,
                               role = "predictor",
                               ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepFindInterval$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_harmonic
#'
#' @description
#'   Add sin and cos terms for harmonic analysis
#'
#' @inheritParams step_scale
#' @param frequency numeric vector - the frequencies of the sin and cos curves
#' @param cycle_size numeric - the period of the sin and cos curves
#' @param starting_value numeric - the starting position of the sin and cos
#'   curves. This may be specified to have more control over the signal phase.
#'
#' @return
#' @export
#'
#' @examples
#' dat <- data.frame(x = 1:10, y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_harmonic(x,
#'                      frequency = 2.0,
#'                      cycle_size = 4.0,
#'                      starting_value = 0.0)
step_harmonic <- function(.rec,
                          terms,
                          frequency = NA_real_,
                          cycle_size = NA_real_,
                          starting_value = 0.0,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepHarmonic$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_intercept
#'
#' @description
#'   Add an intercept term
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#' dat <- data.frame(x = 1:10, y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_intercept()
step_intercept <- function(.rec,
                           terms,
                           role = "predictor",
                           ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepIntercept$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_kernel_filter
#'
#' @description
#'   linearly convolve a kernel with a data series.
#'
#' @param kernel the convolution kernel
#' @param align character center, left or right align the convolution
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
step_kernel_filter <- function(.rec,
                               terms,
                               kernel,
                               align = "center",
                               role = "predictor",
                               ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepIntercept$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_lead_lag
#'
#' @description
#'   Lag or lead a column or columns.  This requires a sorted and regular time
#'   series.
#'
#' @inheritParams step_scale
#' @param lag integer vector - number of samples to lag or lead. Negative
#'   numbers indicate leading a vector.
#' @param n_shift integer - number of values to shift the starting position when
#'   n_subset is not equal to 0. The value of n_shift has to be less than
#'   `n_subset`.
#' @param n_subset integer - spacing between adjacent samples in the result.
#'
#' @return
#' @export
#'
#' @examples
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_lead_lag(x, lag = 1)
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_lead_lag(x, lag = 1, n_subset = 5)
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_lead_lag(x, lag = 1, n_shift = 2, n_subset = 5)
#'
step_lead_lag <- function(.rec,
                          terms,
                          lag,
                          n_shift = 0L,
                          n_subset = 1L,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepLeadLag$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_normalize
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_normalize(x)
#'
step_normalize <- function(.rec,
                           terms,
                           role = "predictor",
                           skip = FALSE,
                           na_rm = TRUE,
                           keep_original_cols = FALSE,
                           ...){

  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepNormalize$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_ols_gap_fill
#'
#' @param recipe Recipe to use for filling gaps
#' @inheritParams step_scale
#'
#' @return
#'
#' @family gap_fill
#'
#' @export
#'
#' @examples
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#'
step_ols_gap_fill <- function(.rec,
                              terms,
                              recipe,
                              role = "predictor",
                              ...){

  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepOlsGapFill$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_ols
#'
#' @description Uses the Eigen C++ library fast versions to generate
#' predictions from different steps.
#'
#'
#' @param recipe Recipe to use for filling gaps
#' @inheritParams step_scale
#'
#' @return
#'
#' @family ols
#'
#' @export
#'
#' @examples
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#'
step_ols <- function(.rec,
                             formula,
                             role = "predictor",
                             do_response = TRUE,
                             do_predict = TRUE,
                             ...){

  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepOls$new,
                        env_list))
}
#' #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' #' @title step_ols_response
#' #'
#' #' @param recipe Recipe to use getting responses from a regression model
#' #' @inheritParams step_scale
#' #'
#' #' @return
#' #'
#' #' @family ols
#' #'
#' #' @export
#' #'
#' #' @examples
#' #' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#' #'
#' #'
#' step_ols_response <- function(.rec,
#'                               formula,
#'                               role = "augment",
#'                               ...){
#'
#'   env_list <- get_function_arguments()
#'   .rec$add_step(do.call(StepOlsResponse$new,
#'                         env_list))
#' }
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_pca
#'
#' @description `StepPca` Does PCA for a set of columns. This currently is an
#' in house function. Use at your own risk!
#'
#' @inheritParams Step
#' @inheritParams recipes::step_pca
#'
#' @return
#' @export
#'
#' @examples
step_pca <- function(.rec,
                     terms,
                     na_rm = TRUE,
                     n_comp = 3,
                     center = TRUE,
                     scale = TRUE,
                     role = "predictor",
                     ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepPca$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_scale
#'
#' @description
#'   Adds a step to scale a data column(s)
#'
#' @param .rec the R6 recipe object.
#' @param terms the unquoted names of the variables to use or a selector
#'   function.  terms replaces the `...` of the recipes package but requires
#'   variables to be included within `c()`.  For example to include variables x
#'   and y you would write `c(x,y)` in the frecipes package.
#' @param role character - the name of the role
#' @param skip logical - should the step be skipped
#' @param na_rm logical - should NA values be removed from calculations
#' @param fun function - the function that is applied to a list or columns of a
#'   data.frame like object.
#' @param n_sd numeric - number of standard deviations for the scaling
#' @param keep_original_cols logical - keep the original columns or replace them
#' @param ... additional arguments
#'
#' @return
#' @export
#'
#' @examples
#'
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_scale(x)
#'
step_scale <- function(.rec,
                       terms,
                       role = "predictor",
                       skip = FALSE,
                       na_rm = TRUE,
                       fun = collapse::fsd,
                       n_sd = 1L,
                       keep_original_cols = FALSE,
                       ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepScale$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_slug_cbp
#'
#' @description
#'   Cooper, Bredehoeft and Papadopulos, 1967 Slug test solution
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#'
#'
step_slug_cbp <- function(.rec,
                          times,
                          radius = 1.0,
                          radius_casing = 1.0,
                          radius_well = 0.15,
                          specific_storage = 1.0e-6,
                          hydraulic_conductivity = 1.0e-4,
                          head_0 = 1.0,
                          thickness = 1.0,
                          n_terms = 16,
                          role = "predictor",
                          ...) {
  times <- substitute(times)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepSlugCbp$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_spline_b
#'
#' @description
#'   generates basis splines.
#'
#' @param internal_knots equivalent to knots from `splines2::bSplines`
#' @param boundary_knots equivalent to Boundary.knots from `splines2::bSplines`
#' @inheritParams splines2::bsp
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#'
#'
step_spline_b <- function(.rec,
                          terms,
                          df = 0L,
                          internal_knots = NULL,
                          boundary_knots = NULL,
                          intercept = FALSE,
                          periodic = FALSE,
                          degree = 3L,
                          role = "predictor",
                          ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepSplineB$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_subset_rows
#'
#' @description
#'   selects rows from output.
#'
#' @param row_numbers integer vector of row numbers to keep.
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#'
#'
step_subset_rows <- function(.rec,
                             terms,
                             row_numbers,
                             role = "modify",
                             ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepSubsetRows$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_transport_fractures_heat
#'
#' @description
#'   Sudicky and Frind 1982 solution adapted for heat. Two parallel fractures.
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
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#'
#'
step_transport_fractures_heat <- function(.rec,
                                          time,
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
  time <- substitute(time)
  distance_fracture <- substitute(distance_fracture)
  distance_matrix <- substitute(distance_matrix)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepTransportFracturesHeat$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_transport_fractures_solute
#'
#' @description
#'   Sudicky and Frind 1982 solution. Two parallel
#' fractures
#'
#' @param time vector elapsed time (t)
#' @param distance_fracture vector distance along fracture (z)
#' @param distance_matrix vector distance into matrix (x)
#' @param concentration_influent vector concentration history (c_in)
#' @param time_influent vector concentration history (t_in)
#' @param concentration_initial double concentration  (c_0)
#' @param fracture_aperture double fracture aperture (2b)
#' @param fracture_spacing double fracture aperture (2B)
#' @param velocity double water velocity in fracture (v)
#' @param dispersivity_longitudinal double longitudinal dispersivity (α_l)
#' @param diffusion double free-water diffusion coefficient (D*)
#' @param sorption_fracture double fracture distribution coefficient (K_f)
#' @param sorption_matrix double matrix distribution coefficient (K_m)
#' @param decay double radioactive half-life for solute (λ)
#' @param density_bulk double dry bulk density (ρ_b)
#' @param porosity double porosity (θ)
#' @param tortuosity double tortuosity (τ)
#' @param n_terms integer number of terms for laplace inversion
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#'
#'
step_transport_fractures_solute <- function(.rec,
                                            time,
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
  time <- substitute(time)
  distance_fracture <- substitute(distance_fracture)
  distance_matrix <- substitute(distance_matrix)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepTransportFracturesSolute$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_transport_ogata_banks
#'
#' @description
#' Ogata, A., Banks, R.B., 1961. A solution of the differential equation of
#' longitudinal dispersion in porous media. U. S. Geol. Surv. Prof. Pap. 411-A.
#' 1-D, infinite source, uniform flow, constant parameters, decay, retardation
#'
#' To have values match the excel sheet
#' https://www.civil.uwaterloo.ca/jrcraig/pdf/OgataBanks.xlsm the decay
#' coefficient needs to be scaled by the retardation coefficient.
#'
#' 1-D
#' infinite source
#' uniform flow
#' constant parameters
#' no decay
#' no retardation
#'
#' @param time vector time
#' @param distance vector x position
#' @param concentration_initial double concentration
#' @param velocity double velocity
#' @param diffusion double diffusion coefficient
#' @param retardation double retardation coefficient
#' @param decay double decay coefficient
#'
#' @inheritParams step_scale
#'
#' @return Ogata-Banks solution for time and distance pairs
#' @export
#'
#' @examples
#'
#'
step_transport_ogata_banks <- function(.rec,
                                       time,
                                       distance,
                                       concentration_initial = 1.0,
                                       velocity = 0.1,
                                       diffusion = 0.1,
                                       retardation = 1.0,
                                       decay = 0.0,
                                       role = "predictor",
                                       ...) {
  time <- substitute(time)
  distance <- substitute(distance)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepTransportOgataBanks$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_vadose_weeks
#'
#' @description
#' Weeks solution
#'
#' @param time vector time
#' @param air_diffusivity double
#' @param thickness double
#' @param precision double
#' @param inverse double
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#'
#'
step_vadose_weeks <- function(.rec,
                              time,
                              air_diffusivity = 0.2,
                              thickness = 40.0,
                              precision = 1e-12,
                              inverse = FALSE,
                              role = "predictor",
                              ...) {
  time <- substitute(time)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepVadoseWeeks$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_varying
#'
#' @description
#'
#' remove columns that only contain a single value.
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#'
#'
step_varying <- function(.rec,
                         terms,
                         role = "predictor",
                         ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepVarying$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# prep -------------------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title prep
#'
#' @description
#'   prep a recipe
#'
#' @inheritParams step_scale
#' @param retain logical - currently not implemented
#'
#' @return
#' @export
#'
#' @examples
#' rec <- recipe(y~x, data = list(x = rnorm(10), y = rnorm(10))) |>
#'        step_scale(x) |>
#'        prep()
prep <- function(.rec, retain = TRUE) {
  .rec$prep(retain)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# bake -------------------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title bake
#'
#' @description
#'   Evaluate the steps and store the recipe results
#'
#' @inheritParams step_scale
#' @inheritParams stats::lm
#' @param type
#'
#' @return
#' @export
#'
#' @examples
#' rec <- recipe(y~x, data = list(x = rnorm(10), y = rnorm(10))) |>
#'        step_scale(x) |>
#'        prep() |>
#'        bake()
bake <- function(.rec, data = NULL) {
  .rec$bake(data = data)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# plate ------------------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title plate
#'
#' @description
#'   Get the results from the recipe. If the recipe hasn't been prepped and
#'   baked, this will do those steps and return the result.
#'
#'
#' @inheritParams step_scale
#' @param type the return type for the recipe (dt = `data.table`, df = `data.frame`,
#' tbl = `tibble`, list = `list`, m = `matrix`)
#'
#' @return
#' @export
#'
#' @examples
#'
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_scale(x) |>
#'        prep() |>
#'        bake() |>
#'        plate()
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_scale(x) |>
#'        plate()
#'
plate <- function(.rec, type = "dt") {
  .rec$plate(type = type)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^





# formula <- as.formula(y~x)
# data <- data.frame(x = as.numeric(1:10000), y = as.numeric(1:10000))
# dat <- data
# frec4 <- frecipes:::recipe(formula, data) |>
#   step_normalize(x) |>
#   prep() |>
#   bake()
#
# bench::mark(
#   rec1 <- recipes::recipe(formula, data) |>
#     recipes::step_scale(x) |>
#     recipes::prep() |>
#     recipes::bake(new_data = NULL),
#   frec2 = Recipe$new(formula = formula, data = data)$
#     add_step(StepScale$new(x))$
#     prep()$
#     bake(),
#   frec1 <- frecipes:::recipe(formula, data) |>
#     step_scale(x) |>
#     prep() |>
#     bake(),
#   frec3 <- frecipes:::recipe(formula, data) |>
#     step_center(x) |>
#     prep() |>
#     bake(),
#   frec4 <- frecipes:::recipe(formula, data) |>
#     step_normalize(x) |>
#     prep() |>
#     bake(),
#   frec5 <- recipes:::recipe(formula, data) |>
#     recipes::step_normalize(x) |>
#     recipes::prep() |>
#     recipes::bake(new_data = NULL),
#   check = FALSE
#   # relative = TRUE
# )


