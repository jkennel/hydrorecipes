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

#' tidal_cooper_1965
#'
#' Cooper Jr, H.H., Bredehoeft, J.D., Papadopulos, I.S. and Bennett, R.R., 1965.
#' The response of well‐aquifer systems to seismic waves. Journal of Geophysical Research, 70(16), pp.3915-3926.
#'
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

  kel   <- kelvin(alpha, nSeq = 1L)

  ker_0 <- kel[['ker_0']]
  kei_0 <- kel[['kei_0']]

  # Equation 28
  e <- 1.0 - (t1 * kei_0) - ((omega)^2 * h_e) / gravity
  f <-        t1 * ker_0

  out <- data.table::data.table(frequency = frequency, period = 1/frequency)
  out[, dimensionless_frequency := transmissivity / (frequency * radius_casing^2)]
  out[, Q := .calc_dimensionless_frequency(omega, height_water, transmissivity/storativity)]
  out[, response := 1.0 / (e + f * 1i )]  # convert to imaginary
  out[, vertical_motion := response * 4.0 * pi^2 * h_e / (period^2 * gravity)]

}

#' tidal_hsieh_1987
#' Solution for estimating transmissivity and storativity from earth tides.
#'
#' Hsieh, P.A., Bredehoeft, J.D. and Farr, J.M., 1987.
#' Determination of aquifer transmissivity from Earth tide analysis. Water resources research, 23(10), pp.1824-1832.
#'
#' @param frequency
#' @param storativity
#' @param transmissivity
#' @param radius_well
#'
#' @return
#' @export
#'
#' @examples
#' data('hsieh_1987_fig_2_3')
#' storativity <- 1e-07
#' transmissivity <- 1e-03
#' radius_well <- 0.05
#' frequency <- 10^seq(-5, 2, by = 0.05)
#' tau   <- 1 / frequency
#' hsieh <- tidal_hsieh_1987(frequency, storativity, transmissivity, radius_well)
#' plot(Mod(response)~dimensionless_frequency, hsieh,
#'  type='l',
#'  log = 'x',
#'  xlim = c(1, 1000))
#' points(response~dimensionless_frequency, hsieh_1987_fig_2_3[variable=='gain' & S == storativity])
#'
#' plot(unwrap(Arg(response)) * 180/pi~dimensionless_frequency, hsieh,
#'  type='l',
#'  log = 'x',
#'  xlim = c(1, 1000),
#'  ylim = c(0, -90))
#' points(response~dimensionless_frequency, hsieh_1987_fig_2_3[variable=='phase' & S == storativity])
#'
tidal_hsieh_1987 <- function(frequency,
                             storativity,
                             transmissivity,
                             radius_well,
                             radius_casing = radius_well) {

  omega <- .calc_omega(frequency)

  # equation 10
  alpha <- .calc_alpha_w(omega, storativity, transmissivity, radius_well)


  kel   <- kelvin(alpha, nSeq = 2)

  ker_0 <- kel[['ker_0']]
  kei_0 <- kel[['kei_0']]

  ker_1 <- kel[['ker_1']]
  kei_1 <- kel[['kei_1']]

  # equations 8 & 9
  denom <- (sqrt(2.0) * alpha * (ker_1^2 + kei_1^2))
  phi <- (-ker_1 - kei_1) / denom
  psi <- (-ker_1 + kei_1) / denom

  # equations 13 & 14
  t1 <- .calc_dimensionless_frequency(omega, radius_casing, transmissivity)
  e  <- 1.0 - t1 * (psi * ker_0 + phi * kei_0)
  f  <-       t1 * (phi * ker_0 - psi * kei_0)

  # return data.table of frequency and response
  out <- data.table::data.table(frequency = rep(frequency, times = length(transmissivity)), psi, phi, e, f)
  out[, dimensionless_frequency := transmissivity / (frequency * radius_casing^2)]
  #out[, a_prime := (1.0 / storativity) * (e^2 + f^2)^(-0.5)]
  out[, response := 1.0 / (e + f * 1i)]  # convert to imaginary

}


#' #' tidal_liu_1989
#' #' Liu, L.B., Roeloffs, E. and Zheng, X.Y., 1989. Seismically induced water level fluctuations in the Wali well, Beijing, China. Journal of Geophysical Research: Solid Earth, 94(B7), pp.9453-9462.
#' #' @param frequency
#' #' @param storativity
#' #' @param transmissivity
#' #' @param thickness_aquifer
#' #' @param height_water
#' #' @param radius_well
#' #' @param radius_casing
#' #' @param gravity
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' #' data('liu_1989_fig_8')
#' #' storativity <- 5e-4
#' #' transmissivity <- 0.5
#' #' radius_well <- 0.117
#' #' tau   <- seq(10, 50, 0.001)
#' #' frequency <- 1/tau
#' #' liu    <- tidal_liu_1989(frequency, storativity, transmissivity, thickness_aquifer = 565, height_water = 92, radius_well, radius_casing = radius_well/2)
#' #' cooper <- tidal_cooper_1965(frequency, storativity, transmissivity, thickness_aquifer = 565, height_water = 92, radius_well, radius_casing = radius_well/2)
#' #' plot(response~period, liu_1989_fig_8, type='p')
#' #' points(Mod(response)~(tau), liu,
#' #'  type='l',
#' #'  ylim = c(0.0, 8), col = 'red')
#' #' points(Mod(response)~tau, cooper,
#' #'  type='l',
#' #'  col = 'blue')
#' #'
#' #' plot(unwrap(Arg(response)) * 180/pi~tau, liu,
#' #'  type='l',
#' #'  log = 'x')
#' #' points(unwrap(Arg(response)) * 180/pi~tau, cooper,
#' #'  type='l',
#' #'  col = 'blue')
#' tidal_liu_1989 <- function(frequency,
#'                            storativity,
#'                            transmissivity,
#'                            thickness_aquifer,
#'                            height_water,
#'                            radius_well,
#'                            radius_casing = radius_well,
#'                            gravity =  9.80665) {
#'
#'   omega <- .calc_omega(frequency)
#'   alpha <- .calc_alpha_w(omega, storativity, transmissivity, radius_well)
#'
#'
#'   kel   <- kelvin(alpha, nSeq = 1)
#'
#'   k_0   <- kel[['k_0']]
#'   # ker_0 <- kel[['ker_0']]
#'   # kei_0 <- kel[['kei_0']]
#'
#'   # Liu equation A13
#'   u <- (thickness_aquifer / transmissivity) * (k_0)
#'
#'   # Liu equation A16
#'   beta <- sqrt((2.0 * 1i * omega) / (radius_well^2 * gravity * u))
#'
#'   # Liu equation A20
#'   beta_term   <- exp(      -beta * thickness_aquifer)
#'   beta_term_2 <- exp(2.0 * -beta * thickness_aquifer)
#'
#'   e <- -(omega^2 / gravity) *
#'     (height_water + (1.0 - beta_term) / (beta * (1.0 + beta_term))) + 1.0
#'
#'   f <- -(1i * omega * u * radius_well^2) *
#'     ((beta * beta_term) / (1 - beta_term_2))
#'
#'
#'   out <- data.table::data.table(frequency)
#'   out[, dimensionless_frequency := transmissivity / (frequency * radius_casing^2)]
#'   out[, Q := .calc_dimensionless_frequency(omega, height_water, transmissivity/storativity)]
#'
#'
#'   # Here we make a modification so that the modulus and phase are correct.
#'   #complex(modulus = Mod(1/(e+f)), argument = Arg(e+f))
#'   out[, response := complex(modulus = Mod(1 / (e + f)),
#'                             argument = -Arg(e - f))]  # convert to imaginary
#'
#' }
#'
#'
#'
#' #' tidal_wang_2018
#' #'
#' #' @param frequency
#' #' @param storativity
#' #' @param transmissivity
#' #' @param k_vertical
#' #' @param thickness_confining
#' #' @param radius_well
#' #' @param radius_casing
#' #'
#' #' @return
#' #' @export
#' #'
#' tidal_wang_2018 <- function(frequency,
#'                             storativity,
#'                             transmissivity,
#'                             k_vertical,
#'                             thickness_confining,
#'                             radius_well,
#'                             radius_casing = radius_well
#' ) {
#'
#'   omega   <- .calc_omega(frequency)
#'   t1      <- 1i * omega * storativity
#'   k_div_b <- k_vertical / thickness_confining
#'   beta    <- sqrt(k_div_b / (transmissivity) + (t1 / transmissivity))
#'   beta_rw <- beta * radius_well
#'   k0      <- Bessel::BesselK(beta_rw, nu = 0, nSeq = 2, expon.scaled = FALSE)
#'
#'   numer <- radius_casing^2 * 1i  * omega * k0[,1]
#'   denom <- radius_well   * 2.0 * transmissivity * beta * k0[,2]
#'
#'   xi    <- 1.0 + (numer / denom)
#'   x0    <- t1 / (xi * (t1 + (k_div_b)))
#'
#'   return(x0)
#'
#' }


