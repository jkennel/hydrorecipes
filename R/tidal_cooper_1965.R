#' kelvin
#' Kelvin functions of the second kind ker and kei and order 0 to 1.
#'
#' @param z value to evaluate the kelvin functions
#'
#' @return data.table of real and imaginary kelvin functions
#'
#' @importFrom data.table ":=" "data.table" "setnames"
#'
#' @export
kelvin <- function(z, nSeq = 2) {

  n  <- 0:(nSeq-1)
  c1 <- exp(1i * pi / 4.0)
  c2 <- matrix(zapsmall(exp(-pi * n * 1i / 2.0), digits = 64),
               ncol = nSeq,
               nrow = length(z),
               byrow = TRUE)

  k  <- collapse::mctl(c2 * Bessel::BesselK(z * c1, nu = 0, nSeq = nSeq))

  nms_k <- c(paste0('k_', n),
             paste0('ker_', n),
             paste0('kei_', n))

  k <- append(k, append(lapply(k, Re),lapply(k, Im)))
  setNames(k, nms_k)

}

# kelvin_2 <- function(z, nSeq = 2) {
#
#   n  <- 0:(nSeq-1)
#   c1 <- exp(1i * pi / 4.0)
#   c2 <- matrix(zapsmall(exp(-pi * n * 1i / 2.0), digits = 64),
#                ncol = nSeq,
#                nrow = length(z),
#                byrow = TRUE)
#
#   k  <- (c2 * Bessel::BesselK(z * c1, nu = 0, nSeq = nSeq))
#
#   collapse::mctl(Re(k))
#
# }
#
# n <- 1e6
# v <- abs(rnorm(n))
# ns <- 2
# bench::mark(kelvin(v, ns), kelvin_2(v,ns), check = FALSE)

#' tidal_cooper_1965
#'
#' @param frequency
#' @param storativity
#' @param transmissivity
#' @param thickness_aquifer
#' @param height_water
#' @param radius_well
#' @param gravity
#'
#' @return
#' @export
#' @examples
#' data('hsieh_1987_fig_2_3')
#' storativity <- 1e-07
#' transmissivity <- 1e-03
#' radius_well <- 0.05
#' frequency <- 10^seq(-5, 2, by = 0.1)
#' tau   <- 1 / frequency
#' cooper <- tidal_cooper_1965(frequency, storativity, transmissivity, thickness_aquifer = 1, height_water = 1, radius_well)
#' plot(Mod(response)~dimensionless_frequency, cooper,
#'  type='l',
#'  log = 'x',
#'  xlim = c(1, 1000))
#' points(response~dimensionless_frequency, hsieh_1987_fig_2_3[variable=='gain' & S == storativity])
#'
#' plot(unwrap(Arg(response)) * 180/pi~dimensionless_frequency, cooper,
#'  type='l',
#'  log = 'x',
#'  xlim = c(1, 1000),
#'  ylim = c(0, -90))
#' points(response~dimensionless_frequency, hsieh_1987_fig_2_3[variable=='phase' & S == storativity])
#'
tidal_cooper_1965 <- function(frequency,
                              storativity,
                              transmissivity,
                              thickness_aquifer,
                              height_water,
                              radius_well,
                              radius_casing = radius_well,
                              gravity =  9.80665) {

  h_e   <- .calc_effective_height(height_water, thickness_aquifer)
  omega <- .calc_omega(frequency)
  alpha <- .calc_alpha_w(omega, storativity, transmissivity, radius_well)
  t1    <- .calc_dimensionless_frequency(omega, radius_casing, transmissivity)

  kel   <- kelvin(alpha, nSeq = 1)

  ker_0 <- kel[['ker_0']]
  kei_0 <- kel[['kei_0']]

  # Equation 28
  e <- 1.0 - (t1 * kei_0) - ((omega)^2 * h_e) / gravity
  f <-        t1 * ker_0

  out <- data.table::data.table(frequency = frequency, period = 1/frequency)
  out[, dimensionless_frequency := transmissivity / (frequency * radius_casing^2)]
  out[, Q := .calc_dimensionless_frequency(omega, height_water, transmissivity/storativity)]
  out[, response := 1.0 / (e + f * 1i )]  # convert to imaginary
  out[, vertical_motion := response * 4 * pi^2 * h_e / (period^2 * gravity)]

}
