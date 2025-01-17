#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#
# recipe -----------------------------------------------------------------------
#
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Create a new R6 recipe. This is analogous to the the list structure that the
#' *recipes* package uses.
#'
#' @description
#' Create a new recipe object.
#' @param formula The model formula. It cannot contain operations.
#' @param data list, data.frame, data.table, tibble of data. They will all
#' be treated as lists.
#'
#' @param ... additional arguments to pass to Recipe$new(). This is currently
#' not used.
#'
#' @return A new R6 `Recipe` object.
#' @export
#'
#' @importFrom collapse fmean fsd fscale fsum fquantile fndistinct flag
#' @importFrom collapse vlengths fnunique
#' @importFrom collapse missing_cases varying rowbind
#' @importFrom collapse qDF qM qF qTBL mctl %!in% "%iin%" "%!iin%"
#' @importFrom collapse pivot
#'
#' @importFrom earthtide calc_earthtide
#' @importFrom R6 R6Class
#' @importFrom Bessel BesselK BesselJ BesselI
#'
#' @importFrom Rcpp sourceCpp
#' @importFrom stats nextn
#' @importFrom stats convolve
#' @importFrom stats spec.pgram
#' @importFrom R6 R6Class
#'
#' @importFrom data.table rleid
#'
#' @useDynLib hydrorecipes, .registration = TRUE
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
#' step_add_noise
#'
#' @description
#'   Add noise.
#'
#' @inheritParams step_scale
#'
#' @param sd the standard deviation of the noise to add
#' @param mean the mean of the noise to add
#' @param fun the random noise generating function
#'
#' @return add noise to a variable
#' @export
#'
#' @examples
#' rec <- recipe(y~x, data = dat) |>
#'        step_add_noise(x) |> plate()
#'
step_add_noise <- function(.rec,
                           terms,
                           sd = 1.0,
                           mean = 0.0,
                           fun = rnorm,
                           role = "predictor",
                           ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepAddNoise$new,
                        modifyList(x = env_list, val = list(...))))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_add_vars
#'
#' @description
#'   Add a variable from the initial (template) data not included in the
#'   recipe creation.
#'
#' @inheritParams step_scale
#'
#' @return recipe with added variables
#' @export
#'
#' @examples
#' dat <- data.frame(x = rnorm(10), y = rnorm(10), z = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_add_vars(z) |>
#'        plate()
#'
step_add_vars <- function(.rec,
                          terms,
                          role = "predictor",
                          ...) {

  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepAddVars$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_aquifer_constant_drawdown
#'
#' @description
#'   Estimates the flows of a constant drawdown test using Jacob-Lohman 1952.
#'
#' @inheritParams step_scale
#' @inheritParams step_aquifer_grf
#' @param radius_well the radius of the well (L)
#' @param n_terms number of terms for laplace solution inversion
#'
#' @return an updated recipe
#'
#' @family aquifer
#'
#' @references
#' Jacob, C.E. and S.W. Lohman, 1952. Nonsteady flow to a well of constant
#'  drawdown in an extensive aquifer, Trans. Am. Geophys. Union, vol. 33,
#'  pp. 559-569.
#'
#' @examples
#' time <- 10^seq(-5, 2, 0.1)
#' form <- formula(time~.)
#' dat <- data.frame(time = time)
#'
#' jl = recipe(formula = form, data = dat) |>
#'   step_aquifer_constant_drawdown(time = time,
#'                                  drawdown = 10,
#'                                  thickness = 10,
#'                                  radius_well = 0.15,
#'                                  specific_storage = 1e-6,
#'                                  hydraulic_conductivity = 1,
#'                                  n_terms = 12L) |>
#'   plate()
#'
#' @export
step_aquifer_constant_drawdown <- function(.rec,
                                           time,
                                           drawdown = 1.0,
                                           thickness = 1.0,
                                           radius_well = 0.15,
                                           specific_storage = 1.0e-6,
                                           hydraulic_conductivity = 1.0e-4,
                                           n_terms = 16L,
                                           role = "predictor",
                                           ...) {
  time <- substitute(time)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAquiferConstantDrawdown$new,
                        modifyList(x = env_list, val = list(...))))

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
#' @inheritParams step_aquifer_constant_drawdown
#' @param time the time for evaluation (t)
#' @param thickness the aquifer thickness (L)
#' @param radius the distance to the observation well (L)
#' @param specific_storage specific storage of aquifer (L/L)
#' @param hydraulic_conductivity the hydraulic conductivity (L/t)
#' @param flow_rate the flow rate from the well (L^3/t)
#' @param flow_dimension aquifer flow dimension (1 = linear, 2 = radial,
#'  3 = spherical)
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
#' time <- 1:2000
#' flow_rate  <- c(rep(0.001, 500),
#'                 rep(0.002, 500),
#'                 rep(0.0, 1000))
#'
#' dat <- data.frame(time, flow_rate)
#'
#' # radial (flow_dimension = 2 Theis)
#' dd_rad <- recipe(time~flow_rate, dat) |>
#'   step_aquifer_grf(time = time,
#'                    flow_rate = flow_rate,
#'                    thickness = 1.0,
#'                    radius = 20,
#'                    specific_storage = 1e-5,
#'                    hydraulic_conductivity = 1e-3,
#'                    flow_dimension = 2) |>
#'   plate()
#'
#'
#' # linear (flow_dimension = 1)
#'
#' dd_lin <- recipe(time~flow_rate, dat) |>
#'   step_aquifer_grf(time = time,
#'                    flow_rate = flow_rate,
#'                    thickness = 1.0,
#'                    radius = 20,
#'                    specific_storage = 1e-5,
#'                    hydraulic_conductivity = 1e-3,
#'                    flow_dimension = 1) |>
#'   plate()
#'
#' # spherical (flow_dimension = 3)
#' dd_sph <- recipe(time~flow_rate, dat) |>
#'   step_aquifer_grf(time = time,
#'                    flow_rate = flow_rate,
#'                    thickness = 1.0,
#'                    radius = 20,
#'                    specific_storage = 1e-5,
#'                    hydraulic_conductivity = 1e-3,
#'                    flow_dimension = 3) |>
#'   plate()
#'
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
                        modifyList(x = env_list, val = list(...))))  # need to add ... to this??


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
#' dat <- data.frame(x = as.numeric(1:100),
#'                   y = rep(0.01, 100))
#' formula <- as.formula(y~x)
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_aquifer_theis(time = x, flow_rate = y) |>
#'   prep() |>
#'   bake()
#'
#' @export
step_aquifer_theis <- function(.rec,
                               time,
                               flow_rate,
                               thickness = 1.0,
                               radius = 100.0,
                               specific_storage = 1.0e-6,
                               hydraulic_conductivity = 1.0e-4,
                               role = "predictor",
                               ...) {
  time <- substitute(time)
  flow_rate <- substitute(flow_rate)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAquiferTheis$new,
                        modifyList(x = env_list, val = list(...))))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_aquifer_theis_aniso
#'
#' @description
#' Generates the drawdown using the Papadopulos 1965 model.
#' Papadopulos, I.S., 1965. Nonsteady flow to a well in an infinite anisotropic
#'   aquifer.
#'
#' @inheritParams step_scale
#' @inheritParams step_aquifer_grf
#' @param distance_x distance in the x direction
#' @param distance_y distance in the y direction
#' @param hydraulic_conductivity_major hydraulic conductivity in the major principal direction
#' @param hydraulic_conductivity_minor hydraulic conductivity in the minor principal direction
#' @param major_axis_angle the orientation of the major principal axis (angle between the x axis and major axis)
#'
#' @return The drawdown using the Papadopulos 1965 model
#'
#' @references
#' Papadopulos, I.S., 1965. Nonsteady flow to a well in an infinite anisotropic
#'   aquifer.
#' Heilweil, V.M. and Hsieh, P.A., 2006. Determining anisotropic transmissivity
#'   using a simplified Papadopulos method. Groundwater, 44(5), pp.749-753.
#'
#' @family aquifer
#'
#' @examples
#' dat <- data.frame(x = as.numeric(1:100),
#'                   y = rep(0.01, 100))
#' formula <- as.formula(y~x)
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_aquifer_theis_aniso(time = x, flow_rate = y) |>
#'   prep() |>
#'   bake()
#'
#' @export
step_aquifer_theis_aniso <- function(.rec,
                                     time,
                                     flow_rate,
                                     thickness = 1.0,
                                     distance_x = 100.0,
                                     distance_y = 100.0,
                                     specific_storage = 1.0e-6,
                                     hydraulic_conductivity_major = 1.0e-4,
                                     hydraulic_conductivity_minor = 1.0e-5,
                                     major_axis_angle = 0.0,
                                     role = "predictor",
                                     ...) {
  time <- substitute(time)
  flow_rate <- substitute(flow_rate)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAquiferTheisAniso$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_aquifer_leaky
#'
#' @description
#' hantush_jacob 1955 solution for a leaky aquifer
#'
#' @inheritParams step_scale
#' @inheritParams step_aquifer_grf
#' @param leakage the leakage defined by hantush (smaller indicates more leaky)
#' @param precision the precision of the solution (default 1e-10)
#'
#' @return The drawdown using the Hantush and Jacob 1955 model
#'
#' @references
#'   Prodanoff, J.H.A., Mansur, W.J. and Mascarenhas, F.C.B., 2006. Numerical
#'   evaluation of Theis and Hantush-Jacob well functions. Journal of hydrology,
#'   318(1-4), pp.173-183. eq: 10, 11, 12
#'
#' Hantush, M.S. and C.E. Jacob, 1955. Non-steady radial flow in an infinite
#'  leaky aquifer, Am. Geophys. Union Trans., vol. 36, no. 1, pp. 95-100.
#'
#' @family aquifer
#'
#' @examples
#' time <- 1:2000
#' flow_rate <- c(rep(0.001, 500),
#'                rep(0.002, 500),
#'                rep(0.0, 1000))
#'
#' # high
#' dat <- data.frame(time = 1:2000, flow_rate = flow_rate)
#' hj_100 <- recipe(flow_rate~time, dat) |>
#'   step_aquifer_leaky(time,
#'                      flow_rate,
#'                      leakage = 100,
#'                      radius = 100,
#'                      storativity = 1e-6,
#'                      transmissivity = 1e-4) |>
#'   plate()
#'
#' # medium
#' hj_200 <- recipe(flow_rate~time, dat) |>
#'   step_aquifer_leaky(time,
#'                      flow_rate,
#'                      leakage = 200,
#'                      radius = 100,
#'                      storativity = 1e-6,
#'                      transmissivity = 1e-4) |>
#'   plate()
#'
#' # low
#' hj_1000 <- recipe(flow_rate~time, dat) |>
#'   step_aquifer_leaky(time,
#'                      flow_rate,
#'                      leakage = 1000,
#'                      radius = 100,
#'                      storativity = 1e-6,
#'                      transmissivity = 1e-4) |>
#'   plate()
#'
#' @export
step_aquifer_leaky <- function(.rec,
                               time,
                               flow_rate,
                               leakage = 100.0,
                               radius = 100.0,
                               storativity = 1e-6,
                               transmissivity = 1e-4,
                               precision = 1e-10,
                               role = "predictor",
                               ...) {
  time <- substitute(time)
  flow_rate <- substitute(flow_rate)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAquiferLeaky$new,
                        modifyList(x = env_list, val = list(...))))

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
#' dat <- data.frame(time = as.numeric(1:100))
#' formula <- as.formula(time~.)
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_aquifer_patch(time = time, flow_rate = 0.01) |>
#'   plate()
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
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_aquifer_patch
#'
#' @description
#' Papadopulos-Cooper 1967 solution for wellbore storage.
#'
#' @inheritParams step_scale
#' @inheritParams step_aquifer_grf
#'
#' @param radius distance from center of well
#' @param radius_casing radius of casing in the interval over which the water
#'   level declines
#' @param radius_well effective radius of well screen or open hole
#'
#'
#' @return The drawdown using the Papadopulos-Cooper model
#'
#' @references
#' Papadopulos, I.S. and H.H. Cooper, 1967. Drawdown in a well of large
#'   diameter, Water Resources Research, vol. 3, no. 1, pp. 241-244.
#'
#'
#' @family aquifer
#'
#' @examples
#' dat <- data.frame(x = 10^seq(-5, 2, length.out = 100),
#'                   y = rep(0.01, 100))
#' formula <- as.formula(y~x)
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_aquifer_wellbore_storage(time = x,
#'                                 flow_rate = 10.0,
#'                                 hydraulic_conductivity = 10.0,
#'                                 specific_storage = 1e-4) |>
#'   prep() |>
#'   bake()
#'
#' @export
step_aquifer_wellbore_storage <- function(.rec,
                                          time,
                                          flow_rate = 1.0,
                                          radius = 0.15,
                                          radius_casing = 0.15,
                                          radius_well = 0.15,
                                          thickness = 1.0,
                                          specific_storage = 1.0e-6,
                                          hydraulic_conductivity = 1.0e-4,
                                          n_terms = 12L,
                                          role = "predictor",
                                          ...) {
  time <- substitute(time)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepAquiferWellboreStorage$new,
                        modifyList(x = env_list, val = list(...))))

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
#' data(kennel_2020)
#'
#' clarks <- recipe(wl~., kennel_2020) |>
#'   step_baro_clark(wl, baro, lag_space = 1) |> # 1 minutes (every minute differences)
#'   step_baro_clark(wl, baro, lag_space = 60) |> # 60 minutes (hourly differences)
#'   step_baro_clark(wl, baro, lag_space = 1440) |> # 1440 minutes (daily differences)
#'   prep() |>
#'   bake()
#'
#' clarks$get_step_data("barometric_efficiency")
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
                        modifyList(x = env_list, val = list(...))))

}



# Hussein
# Determination of Fluid Flow Properties From the Response of Water Levels in Wells to Atmospheric Loading
#' step_baro_frequency_semi_confined
#'
#' @description
#' Work in Progress: Do not use. Rojstaczer 1988 solution
#'
#' step_baro_frequency_semi_confined
#'
#' @inheritParams step_scale
#'
#' @param frequency
#' @param radius_well
#' @param transmissivity
#' @param storage_confining
#' @param storage_aquifer
#' @param diffusivity_confining
#' @param diffusivity_vadose
#' @param thickness_confining
#' @param thickness_vadose
#' @param loading_efficiency
#' @param attenuation
#'
#' @return complex response vector in frequency domain
#'
#' @export
step_baro_frequency_semi_confined <- function(.rec,
                                              frequency,
                                              radius_well,
                                              transmissivity,
                                              storage_confining,
                                              storage_aquifer,
                                              diffusivity_confining,
                                              diffusivity_vadose,
                                              thickness_confining,
                                              thickness_vadose,
                                              loading_efficiency,
                                              attenuation,
                                              role = "predictor",
                                              ...) {
  frequency <- substitute(frequency)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepBaroFrequencySemiConfined$new,
                        modifyList(x = env_list,
                                   val = list(...))))

}



#' step_baro_frequency_semi_confined
#'
#' @description
#' Work in Progress: Do not use. Rojstaczer 1988 solution
#'
#' step_baro_frequency_unconfined
#'
#' @inheritParams step_scale
#' @inheritParams step_baro_frequency_semi_confined
#'
#' @param specific_yield
#' @param k_vertical
#' @param diffusivity_vertical
#' @param thickness_saturated_well
#' @param thickness_vadose
#' @param thickness_aquifer
#'
#' @return complex response vector in frequency domain
#'
#' @export
step_baro_frequency_unconfined <- function(.rec,
                                           frequency,
                                           radius_well,
                                           storage_aquifer,
                                           specific_yield,
                                           k_vertical,
                                           diffusivity_vertical,
                                           diffusivity_vadose,
                                           thickness_saturated_well,
                                           thickness_vadose,
                                           thickness_aquifer,
                                           loading_efficiency,
                                           attenuation,
                                           role = "predictor",
                                           ...) {
  frequency <- substitute(frequency)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepBaroFrequencyUnconfined$new,
                        modifyList(x = env_list,
                                   val = list(...))))

}

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_baro_clark
#'
#' @description
#' Least squares solution for calculating barometric efficiency
#'
#' @inheritParams step_baro_clark
#'
#' @param differences \code{numeric vector} number of samples between differences
#'
#' @return barometric efficiency using least squares
#'
#' @family barometric
#'
#'
#' @examples
#'
#' data(kennel_2020)
#'
#' least_squares <- recipe(wl~., kennel_2020) |>
#'   step_baro_least_squares(wl, baro) |> # 1 minutes (every minute differences)
#'   step_baro_least_squares(wl, baro, lag_space = 1440, differences = TRUE) |> # 1440 minutes (daily differences)
#'   prep() |>
#'   bake()
#'
#' least_squares$get_step_data("barometric_efficiency")
#'
#' @export
step_baro_least_squares <- function(.rec,
                                    water_level,
                                    barometric_pressure,
                                    lag_space = 1L,
                                    inverse = FALSE,
                                    differences = FALSE,
                                    role = "augment",
                                    ...) {
  water_level <- substitute(water_level)
  barometric_pressure <- substitute(barometric_pressure)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepBaroLeastSquares$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' step_baro_harmonic
#'
#' @description
#' Estimate the static barometric efficiency using harmonic methods
#'   (Ratio, Acworth, Rau, Transfer)
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
#' @return \code{double} barometric efficiency using different methods
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
#' data("kennel_2020")
#' library(data.table)
#' library(collapse)
#'
#' formula <- as.formula(wl~.)
#' frec = recipe(formula = formula, data = kennel_2020) |>
#'  step_baro_harmonic(datetime,
#'                     wl,
#'                     baro,
#'                     et,
#'                     inverse = FALSE)
#' @export
step_baro_harmonic <- function(.rec,
                               time,
                               water_level,
                               barometric_pressure,
                               earth_tide,
                               frequency = c(1.9324, 2.0),
                               cycle_size = 86400,
                               start = 0.0,
                               inverse = FALSE,
                               role = "augment",
                               ...) {
  time <- substitute(time)
  water_level <- substitute(water_level)
  barometric_pressure <- substitute(barometric_pressure)
  earth_tide <- substitute(earth_tide)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepBaroHarmonic$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_center
#' @description
#'   Adds a step to center a data column(s)
#'
#' @inheritParams step_scale
#'
#' @return an updated recipe
#'
#' @export
#'
#' @examples
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_center(x) |>
#'        prep() |>
#'        bake()
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
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_check_na
#'
#' @description
#'   Check columns for NA
#'
#' @inheritParams step_scale
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_check_na(x) |>
#'        prep() |>
#'        bake()
#'
step_check_na <- function(.rec,
                          terms,
                          role = "check",
                          ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepCheckNA$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_check_spacing
#'
#' @description
#'   Check the spacing of a variable
#'
#' @inheritParams step_scale
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_check_spacing(x) |>
#'        prep() |>
#'        bake()
#'
step_check_spacing <- function(.rec,
                               terms,
                               role = "check",
                               ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepCheckSpacing$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_check_spacing
#'
#' @description
#'   Check the spacing of a variable
#'
#' @inheritParams step_scale
#'
#' @param data \code{variable} unquoted data column name
#' @param compare \code{variable} unquoted column name for comparison
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' data("kennel_2020")
#'
#' kennel_2020[1e4, wl := 13.36]
#'
#' frec = recipe(wl~baro, data = kennel_2020) |>
#'  step_compare_columns(data = wl, compare = baro, n_sd = 15) |>
#'  prep() |>
#'  bake()
#'
step_compare_columns <- function(.rec,
                                 data,
                                 compare,
                                 role = "add",
                                 n_sd = 4,
                                 na_rm = TRUE,
                                 ...) {
  data <- substitute(data)
  compare <- substitute(compare)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepCompareColumns$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_convolve_gamma
#'
#' @description
#'   linearly convolve a gamma kernel with a data series.
#'
#' @param amplitude amplitude
#' @param k shape
#' @param theta scale
#'
#' @inheritParams step_kernel_filter
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' formula <- as.formula(x~y+z)
#' rows <- 1e4
#'
#' dat <- data.frame(x = rep(1, rows),
#'                   y = 1:rows,
#'                   z = cumsum(rnorm(rows)))
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_convolve_gamma(z, amplitude = 1, theta = 1, k = 1) |>
#'   plate("tbl")
#'
step_convolve_gamma <- function(.rec,
                                terms,
                                amplitude,
                                k,
                                theta,
                                align = "right",
                                max_length = Inf,
                                role = "predictor",
                                ...) {

  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepConvolveGamma$new,
                        modifyList(x = env_list, val = list(...))))


}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_convolve_exponential
#'
#' @description
#'   linearly convolve a gamma kernel with a data series.
#'
#' @param amplitude amplitude
#' @param theta scale
#'
#' @inheritParams step_kernel_filter
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' formula <- as.formula(x~y+z)
#' rows <- 1e4
#'
#' dat <- data.frame(x = rep(1, rows),
#'                   y = 1:rows,
#'                   z = cumsum(rnorm(rows)))
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_convolve_gamma(z, amplitude = 1, theta = 1, k = 1) |>
#'   plate("tbl")
#'
step_convolve_exponential <- function(.rec,
                                      terms,
                                      amplitude,
                                      theta,
                                      align = "right",
                                      max_length = Inf,
                                      role = "predictor",
                                      ...) {

  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepConvolveExponential$new,
                        modifyList(x = env_list, val = list(...))))


}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_cross_correlation
#'
#' @description
#'   Calculate the autocorrelation function or cross-correlation
#'
#' @inheritParams step_scale
#'
#' @references
#' @return an updated recipe
#' @export
#'
#' @examples
#' formula <- as.formula(y~x)
#' rows <- 1e4
#'
#' dat <- data.frame(x = rnorm(rows),
#'                   y = as.numeric(1:rows))
#'
#' frec = recipe(formula = formula, data = dat) |>
#'  step_cross_correlation(x)
#'
#' frec = recipe(formula = formula, data = dat) |>
#'  step_cross_correlation(c(x, y))
#'
step_cross_correlation <- function(.rec,
                                   terms,
                                   lag_max = 100,
                                   role = "predictor",
                                   ...) {

  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepCrossCorrelation$new,
                        modifyList(x = env_list, val = list(...))))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_distributed_lag
#'
#' @description
#'   Generates distributed lag vectors. For regular spaced lags this uses an FFT
#'   based method which is faster and more memory efficient.
#'
#' @inheritParams step_scale
#'
#' @references
#' Gasparrini, A., 2011. Distributed Lag Linear and Non-Linear Models in R:
#'   The Package dlnm. Journal of statistical software 465 43, 1–20.
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#' formula <- as.formula(y~x)
#' rows <- 1e4
#'
#' dat <- data.frame(x = rnorm(rows),
#'                   y = as.numeric(1:rows),
#'                   z = rnorm(rows))
#'
#' frec = recipe(formula = formula, data = dat) |>
#'  step_distributed_lag(x, knots = hydrorecipes:::log_lags_arma(6, 800))
#'
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
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_drop_columns
#'
#' @description
#'   Removes regressors from recipe.
#'
#' @inheritParams step_scale
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#' rows <- 20
#' formula <- as.formula(y~x)
#'
#' dat <- data.frame(x = rnorm(rows),
#'                   y = as.numeric(1:rows),
#'                   z = rnorm(rows))
#'
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_drop_columns(x)
step_drop_columns <- function(.rec,
                              terms,
                              role = "modify",
                              ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepDropColumns$new,
                        modifyList(x = env_list, val = list(...))))

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
#' @return an updated recipe
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
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_earthtide
#'
#' @description
#'   Generate synthetic Earth tide waves and wave groups.
#'
#' @inheritParams step_scale
#' @inheritParams earthtide::calc_earthtide
#' @param interp_factor calculate for every \code{interp_factor} samples. This
#'   may be faster while being more accurate than adjusting the cutoff.
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#' data(kennel_2020)
#' latitude     <- 34.23411                           # latitude
#' longitude    <- -118.678                           # longitude
#' elevation    <- 500                                # elevation
#' cutoff       <- 1e-5                               # cutoff
#' catalog      <- 'ksm04'                            # hartmann wenzel catalog
#' astro_update <- 300                                # how often to update astro parameters
#' method       <- 'volume_strain'                    # which potential to calculate
#'
#' wave_groups_dl <- as.data.table(earthtide::eterna_wavegroups)
#' wave_groups_dl <- na.omit(wave_groups_dl[time == '1 month'])
#' wave_groups_dl <- wave_groups_dl[wave_groups_dl$start > 0.5,]
#' wave_groups_dl <- wave_groups_dl[, list(start, end)]
#' ngr <- nrow(wave_groups_dl)
#'
#' rec <- recipe(wl~baro+datetime, data = kennel_2020) |>
#'   step_earthtide(datetime,
#'                  wave_groups = wave_groups_dl,
#'                  latitude = latitude,
#'                  longitude = longitude,
#'                  elevation = elevation,
#'                  cutoff = cutoff,
#'                  catalog = catalog)
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
                           astro_update = 1L,
                           interp_factor = 1L,
                           role = "predictor",
                           ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepEarthtide$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_fft_coherence
#'
#' @description
#'   estimates the coherence between terms.
#'
#' @inheritParams step_scale
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' dat <- data.frame(x = rnorm(200),
#'                   y = rnorm(200),
#'                   z = rnorm(200))
#'
#' formula <- as.formula(.~x+y)
#' frec = recipe(formula = formula, data = dat) |>
#'   step_fft_pgram(c(x, y,z)) |>
#'   step_fft_coherence() |>
#'   plate("df")
#'
step_fft_coherence <- function(.rec,
                               terms,
                               role = "augment",
                               ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepCoherence$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_fft_pgram
#'
#' @description
#'   Periodgrams and cross-periodograms using a method similar to
#'   \code{stats::spec.pgram}.
#'
#' @inheritParams step_scale
#' @inheritParams stats::spec.pgram
#' @param lst \code{logical} return a list?
#' @param pad_fft \code{logical} Zero pad the list for faster FFT calculation?
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#' formula <- as.formula(y~.)
#'
#' dat <- data.frame(x = rnorm(200),
#'                   y = rnorm(200))
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_fft_pgram(c(x,y))
#'
step_fft_pgram <- function(.rec,
                           terms,
                           spans = 3,
                           detrend = TRUE,
                           demean = TRUE,
                           lst = TRUE,
                           taper = 0.1,
                           pad_fft = TRUE,
                           time_step = 1,
                           role = "predictor",
                           ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepPgram$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_fft_welch
#'
#' @description
#'  calculates the periodogram (estimate of spectral density) using
#'  Welch's method.
#'
#' @inheritParams step_scale
#' @param length_subset length of fft section
#' @param overlap amount of overlap
#' @param window window weights
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#' formula <- as.formula(y~.)
#'
#' dat <- data.frame(x = rnorm(200),
#'                   y = rnorm(200))
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_fft_welch(c(x,y), length_subset = 10, window = window_rectangle(10))
step_fft_welch <- function(.rec,
                           terms,
                           length_subset,
                           overlap = 0.8,
                           window,
                           time_step = 1.0,
                           role = "augment",
                           ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepWelch$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_fft_transfer_pgram
#'
#' @description
#'  Calculates the transfer function with pgram results.
#'
#' @inheritParams step_scale
#' @inheritParams stats::spec.pgram
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#' data(kennel_2020)
#'
#' form <- as.formula("wl~.")
#'
#' rec <- recipe(form, kennel_2020) |>
#'        step_fft_transfer_pgram(c(wl, baro, et), spans = 3) |>
#'        plate()
#'
step_fft_transfer_pgram <- function(.rec,
                                    terms,
                                    spans = 3,
                                    detrend = TRUE,
                                    demean = TRUE,
                                    taper = 0.1,
                                    time_step = 1.0,
                                    formula = NULL,
                                    role = "augment",
                                    ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepTransferPgram$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_fft_transfer_experimental
#'
#' @description
#'  Calculates the transfer function with pgram results.
#'
#' @inheritParams step_scale
#' @inheritParams stats::spec.pgram
#' @param power spacing for the groups
#' @param n_groups number of results
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#' data(kennel_2020)
#'
#' form <- as.formula("wl~.")
#'
#' rec <- recipe(form, kennel_2020) |>
#'        step_fft_transfer_experimental(c(wl, baro, et), spans = 3) |>
#'        plate()
#'
step_fft_transfer_experimental <- function(.rec,
                                           terms,
                                           spans = 3,
                                           detrend = TRUE,
                                           demean = TRUE,
                                           taper = 0.1,
                                           # power = 3,
                                           n_groups = 200,
                                           time_step = 1.0,
                                           formula = NULL,
                                           role = "augment",
                                           ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepTransferExperimental$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_fft_transfer_welch
#'
#' @description
#'  calculates the transfer function with Welch's results.
#'
#' @inheritParams step_scale
#' @inheritParams step_fft_welch
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#' data(kennel_2020)
#' form <- as.formula("wl~.")
#'
#'   rec <- recipe(form, kennel_2020) |>
#'   step_fft_transfer_welch(c(wl, baro, et),
#'                           length_subset = 1440*8 + 1,
#'                           overlap = 0.6,
#'                           window = window_nuttall(1440*8+1)) |>
#'   plate()
#'
step_fft_transfer_welch <- function(.rec,
                                    terms,
                                    length_subset,
                                    overlap = 0.8,
                                    window,
                                    time_step = 1.0,
                                    formula = NULL,
                                    role = "augment",
                                    ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepTransferWelch$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_find_interval
#'
#' @description
#' divides a series into intervals and then performs dummy encoding.
#'
#' @param vec a vector of break points
#' @inheritParams step_scale
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' formula <- as.formula(y~x)
#'
#' rows = 200
#' dat <- data.frame(x = rnorm(rows),
#'                   y = 1:rows,
#'                   z = rnorm(rows))
#'
#' frec1 = recipe(formula = formula, data = dat) |>
#'   step_find_interval(x, vec = c(-0.1, 0.0, 0.1)) |>
#'   plate("tbl")
#'
step_find_interval <- function(.rec,
                               terms,
                               vec,
                               role = "augment",
                               ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepFindInterval$new,
                        modifyList(x = env_list, val = list(...))))

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
#' @return an updated recipe
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
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_intercept
#'
#' @description
#'   Add an intercept term
#'
#' @inheritParams step_scale
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#' dat <- data.frame(x = 1:10, y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_intercept()
step_intercept <- function(.rec,
                           terms,
                           value = 1.0,
                           role = "predictor",
                           ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepIntercept$new,
                        modifyList(x = env_list, val = list(...))))

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
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' formula <- as.formula(x~y+z)
#' rows <- 1e4
#'
#' dat <- data.frame(x = rep(1, rows),
#'                   y = 1:rows,
#'                   z = cumsum(rnorm(rows)))
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_kernel_filter(z, kernel = list(rep(1, 1001)/1001), align = "center") |>
#'   plate("tbl")
#'
step_kernel_filter <- function(.rec,
                               terms,
                               kernel,
                               align = "center",
                               role = "predictor",
                               ...) {

  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepKernelFilter$new,
                        modifyList(x = env_list, val = list(...))))


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
#' @return an updated recipe
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
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_multiply
#'
#' @inheritParams step_scale
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#' dat <- data.frame(x = rnorm(10), y = rnorm(10))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_multiply(x, value = 4)
#'
step_multiply <- function(.rec,
                          terms,
                          values = 1.0,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...){

  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepMultiply$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_normalize
#'
#' @inheritParams step_scale
#'
#' @return an updated recipe
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
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_ols_gap_fill
#'
#' @param recipe Recipe to use for filling gaps
#' @inheritParams step_scale
#'
#' @return an updated recipe
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
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_ols
#'
#' @description Uses the Eigen C++ library fast versions to generate
#' predictions and coefficients from a recipe.
#'
#'
#' @inheritParams step_scale
#' @param do_response \code{logical} calculate and return the responses?
#' @param formula formula for the regression
#'
#' @return an updated recipe
#'
#' @family ols
#'
#' @export
#'
#' @examples
#' data("kennel_2020")
#' kennel_2020[, datetime := as.numeric(datetime)]
#' formula <- as.formula(wl~.)
#' n_knots <- 12
#' deg_free <- 27
#' max_lag <- 1 + 720
#'
#' frec = recipe(formula = formula, data = unclass(kennel_2020)) |>
#'   step_distributed_lag(baro, knots = hydrorecipes:::log_lags_arma(n_knots, max_lag)) |>
#'   step_spline_b(datetime, df = deg_free, intercept = FALSE) |>
#'   step_intercept() |>
#'   step_drop_columns(baro) |>
#'   step_drop_columns(datetime) |>
#'   step_ols(formula) |>
#'   prep() |>
#'   bake()
step_ols <- function(.rec,
                     formula,
                     role = "predictor",
                     do_response = TRUE,
                     # do_predict = TRUE,
                     ...){

  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepOls$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_nls
#'
#' @description Uses the Eigen C++ library fast versions to generate
#' predictions and coefficients from a recipe.
#'
#'
#' @inheritParams step_scale
#' @param do_response \code{logical} calculate and return the responses?
#' @param do_predict calculate and return the predictions?
#' @param formula formula for the regression
#'
#' @return an updated recipe
#'
#' @family ols
#'
#' @export
#'
#' @examples
#' data("kennel_2020")
#' kennel_2020[, datetime := as.numeric(datetime)]
#' formula <- as.formula(wl~.)
#' n_knots <- 12
#' deg_free <- 27
#' max_lag <- 1 + 720
#'
#' frec = recipe(formula = formula, data = unclass(kennel_2020)) |>
#'   step_distributed_lag(baro, knots = hydrorecipes:::log_lags_arma(n_knots, max_lag)) |>
#'   step_spline_b(datetime, df = deg_free, intercept = FALSE) |>
#'   step_intercept() |>
#'   step_drop_columns(baro) |>
#'   step_drop_columns(datetime) |>
#'   step_ols(formula) |>
#'   prep() |>
#'   bake()
step_nls <- function(.rec,
                     formula,
                     role = "predictor",
                     algorithm = "lm",
                     n_subset = 1L,
                     n_shift = 0L,
                     # do_response = TRUE,
                     # do_predict = TRUE,
                     ...){

  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepNls$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_pca
#'
#' @description `StepPca` Does PCA for a set of columns. This currently is an
#' in house function. Use at your own risk!
#'
#' @inheritParams step_scale
#' @inheritParams recipes::step_pca
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' set.seed(1)
#'
#' formula <- as.formula(x~a+b+d+e+f+g)
#' rows <- 1000
#'
#' dat <- data.frame(x = rnorm(rows),
#'                   a = rnorm(rows),
#'                   b = rnorm(rows),
#'                   d = rnorm(rows),
#'                   e = rnorm(rows),
#'                   f = rnorm(rows),
#'                   g = rnorm(rows))
#'
#' rec  = recipe(formula = formula, data = dat) |>
#'   step_pca(all_numeric()) |>
#'   plate()
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
                        modifyList(x = env_list, val = list(...))))

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
#'   and y you would write `c(x,y)` in the hydrorecipes package.
#' @param role character - the name of the role
#' @param skip logical - should the step be skipped
#' @param na_rm logical - should NA values be removed from calculations
#' @param fun function - the function that is applied to a list or columns of a
#'   data.frame like object.
#' @param n_sd numeric - number of standard deviations for the scaling
#' @param keep_original_cols logical - keep the original columns or replace them
#' @param ... additional arguments
#'
#' @return an updated recipe
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
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_slug_cbp
#'
#' @description
#'   Cooper, Bredehoeft and Papadopulos, 1967 Slug test solution
#'
#' @inheritParams step_scale
#' @inheritParams step_aquifer_constant_drawdown
#'
#' @param radius distance from line source or center of well
#' @param radius_casing radius of casing
#' @param radius_well radius of well screen
#' @param head_0 initial displacement
#'
#' @return an updated recipe
#'
#' @family slug
#'
#' @references
#' Cooper, H.H., J.D. Bredehoeft and S.S. Papadopulos, 1967. Response of a
#'  finite-diameter well to an instantaneous charge of water, Water Resources
#'  Research, vol. 3, no. 1, pp. 263-269.
#'
#' @export
#'
#' @examples
#' # check vs. CBP 1967 table 1
#' times    <- rep(c(1.0, 2.15, 4.64), 4) * 10^(rep(c(-3, -2, -1, 0), each = 3))
#' dat <- list(x = times)
#'
#' formula = formula(x~.)
#'
#' frec1 = recipe(formula = formula, data = dat) |>
#'   step_slug_cbp(
#'     times = x,
#'     radius = 1.0,
#'     radius_casing = 1.0,
#'     radius_well = 1.0,
#'     specific_storage = 1e-1,
#'     hydraulic_conductivity = 1.0,
#'     thickness = 1.0,
#'     head_0 = 1.0,
#'     n_terms = 12L
#'   ) |>
#'   plate("dt")
step_slug_cbp <- function(.rec,
                          times,
                          radius = 1.0,
                          radius_casing = 0.15,
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
                        modifyList(x = env_list, val = list(...))))

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
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' formula <- as.formula(x~y+z)
#' rows <- 1e5
#'
#' dat <- data.frame(x = rnorm(rows),
#'                   y = 1:rows,
#'                   z = cumsum(rnorm(rows)))
#' ik <- collapse::fquantile(dat$x, probs = seq(0, 1, 0.1))
#' bk <- ik[c(1, length(ik))]
#' ik <- ik[-c(1, length(ik))]
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_spline_b(x, df = 11L, intercept = FALSE)  |>
#'  plate("tbl")
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
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_spline_n
#'
#' @description
#'   generates basis splines.
#'
#' @param internal_knots equivalent to knots from `splines2::bSplines`
#' @param boundary_knots equivalent to Boundary.knots from `splines2::bSplines`
#' @inheritParams splines2::nsp
#' @inheritParams step_scale
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' formula <- as.formula(x~y+z)
#' rows <- 1e5
#'
#' dat <- data.frame(x = rnorm(rows),
#'                   y = 1:rows,
#'                   z = cumsum(rnorm(rows)))
#' ik <- collapse::fquantile(dat$x, probs = seq(0, 1, 0.1))
#' bk <- ik[c(1, length(ik))]
#' ik <- ik[-c(1, length(ik))]
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_spline_n(x, df = 11L, intercept = FALSE)  |>
#'  plate("tbl")
#'
step_spline_n <- function(.rec,
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
  .rec$add_step(do.call(StepSplineN$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_subset_na_omit
#'
#' @description
#'   selects rows from output.
#'
#' @inheritParams step_scale
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
step_subset_na_omit <- function(.rec,
                             terms,
                             role = "modify",
                             ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepSubsetNAOmit$new,
                        modifyList(x = env_list, val = list(...))))

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
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' dat <- data.frame(x = as.numeric(1:200),
#' y = rnorm(200))
#' formula <- as.formula(y~x)
#'
#' frec1 = recipe(formula = formula, data = dat) |>
#'   step_subset_rows(y, row_numbers = c(1, 5, 10)) |>
#'   plate("dt")
#'
step_subset_rows <- function(.rec,
                             terms,
                             row_numbers,
                             role = "modify",
                             ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepSubsetRows$new,
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_subset_sample
#'
#' @description
#'   selects rows from output.
#'
#' @param size number of samples.
#' @inheritParams step_scale
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#'
step_subset_sample <- function(.rec,
                             terms,
                             size,
                             role = "modify",
                             ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepSubsetSample$new,
                        modifyList(x = env_list, val = list(...))))
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
#' @family transport
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' formula <- as.formula(~time+z+x)
#'
#' dat <- setDT(expand.grid(10^(3:8),
#'                          seq(0.0, 100, 1),
#'                          c(0.0, 0.05)))
#'
#' names(dat) <- c("time", "z", "x")
#'
#'
#' frec1 = recipe(formula = formula, data = dat) |>
#'   step_transport_fractures_heat(time = time,
#'                                 distance_fracture = z,
#'                                 distance_matrix = x) |>
#'   plate("dt")
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
                        modifyList(x = env_list, val = list(...))))

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
#' @family transport
#'
#' @references
#' Sudicky, E.A., Frind, E.O., Contaminant transport in fractured porous media:
#'  Analytical solutions for a system of parallel fractures, December 1982
#'  https://doi.org/10.1029/WR018i006p01634
#'
#' @return an updated recipe
#' @export
#'
#' @examples
#' formula <- as.formula(~time+z+x)
#'
#' dat <- setDT(expand.grid(10^(3:8),
#'                          seq(0.0, 10, 1),
#'                          c(0.0)))
#'
#' names(dat) <- c("time", "z", "x")
#'
#' frec1 = recipe(formula = formula, data = dat) |>
#'   step_transport_fractures_solute(time = time,
#'                                   distance_fracture = z,
#'                                   distance_matrix = x) |>
#'   plate("dt")
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
                        modifyList(x = env_list, val = list(...))))

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
#' Care must be taken so that input values do not lead to NaN. -Need to fix this.
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
#' @family transport
#'
#' @references
#' Ogata, A., Banks, R.B., 1961. A solution of the differential equation of
#'   longitudinal dispersion in porous media. U. S. Geol. Surv. Prof. Pap. 411-A.
#'   1-D, infinite source, uniform flow, constant parameters, decay, retardation
#'
#' @return Ogata-Banks solution for time and distance pairs
#' @export
#'
#' @examples
#'
#' formula <- as.formula(y~x)
#' rows <- 100
#'
#' dat <- data.frame(expand.grid(as.numeric(1:rows), as.numeric(1:10)))
#' names(dat) <- c('x', 'y')
#' frec1 = recipe(formula = formula, data = dat) |>
#'   step_transport_ogata_banks(time = x, distance = y) |>
#'   plate("dt")
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
                        modifyList(x = env_list, val = list(...))))

}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_vadose_weeks
#'
#' @description
#' Weeks 1979 solution
#'
#' @param time vector time elapsed time from start of pressure change
#' @param air_diffusivity double vadose zone air diffusivity
#' @param thickness double vadose zone thickness
#' @param precision double stop the sum when precision is reached
#' @param inverse double whether the response is invers
#'
#' @inheritParams step_scale
#'
#' @return an updated recipe
#' @family vadose
#'
#' @references
#' Weeks, E.P., 1979. Barometric fluctuations in wells tapping deep unconfined
#'   aquifers. Water Resources Research, 15(5), pp.1167-1176.
#'
#' @export
#'
#'
#' @examples
#' formula <- as.formula(y~x)
#'
#' n <- 100
#' dat <- data.frame(x = as.numeric(1:rows),
#'                   y = as.numeric(1:rows))
#'
#' frec1 = recipe(formula = formula, data = dat) |>
#'   step_vadose_weeks(time = x,
#'                     air_diffusivity = 0.8,
#'                     thickness = 5,
#'                     precision = 1e-12) |>
#'   plate("dt")
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
                        modifyList(x = env_list, val = list(...))))

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
#' @return an updated recipe
#' @export
#'
#' @examples
#'
#' formula <- as.formula(y~x+z)
#' rows <- 1000
#' dat <- data.frame(x = rep(1, rows),
#'                   y = 1:rows,
#'                   z = rnorm(rows))
#'
#' frec = recipe(formula = formula, data = dat) |>
#'   step_varying(c(x, y, z)) |>
#'   plate("tbl")
step_varying <- function(.rec,
                         terms,
                         role = "predictor",
                         ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepVarying$new,
                        modifyList(x = env_list, val = list(...))))

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
#' @return an updated recipe
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
#' @return an updated recipe
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
#' @return an updated recipe
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
plate <- function(.rec, type = "dt", ...) {
  .rec$plate(type = type, ...)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^





# formula <- as.formula(y~x)
# data <- data.frame(x = as.numeric(1:10000), y = as.numeric(1:10000))
# dat <- data
# frec4 <- hydrorecipes:::recipe(formula, data) |>
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
#   frec1 <- hydrorecipes:::recipe(formula, data) |>
#     step_scale(x) |>
#     prep() |>
#     bake(),
#   frec3 <- hydrorecipes:::recipe(formula, data) |>
#     step_center(x) |>
#     prep() |>
#     bake(),
#   frec4 <- hydrorecipes:::recipe(formula, data) |>
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


