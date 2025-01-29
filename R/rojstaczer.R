# THESE ARE WORKS IN PROGRESS DO NOT USE

# helpers -----------------------------------------------------------------
#' unwrap
#' Removes the large phase shifts that can be associated with using Arg or atan2.
#'
#' @param phase the phase in radians
#'
#' @return unwrapped vector
#' @export
#'
unwrap <- function(phase) {

  d <- c(0, -diff(phase))
  s <- 2 * pi * (((d > pi) > 0) - ((d < -pi) > 0))

  phase + cumsum(s)

}

# .calc_omega -------------------------------------------------------------
# convert from frequency (cycles per time) to angular frequency (radians per time)
.calc_omega <- function(frequency) {
  2 * pi * frequency
}


# .calc_W -----------------------------------------------------------------
# Rojstaczer, 1988 Determination of Fluid Flow Properties From the Response of
# Water Levels in Wells to Atmospheric Loading
# Equation 14
.calc_W <- function(omega, radius, transmissivity) {
  (radius^2 * omega) / transmissivity
}


# .calc_alpha_w -----------------------------------------------------------
# Cooper, Bredehoeft, Papadopulos, Bennett, 1965: Equation 19
# Hsieh, 1987: Equation 10
# Liu, Roeloffs, Zheng, 1989: Equation 4
.calc_alpha_w <- function(omega, storativity, transmissivity, radius_well) {
  sqrt((omega * storativity) / transmissivity) * radius_well
}


# .calc_dimensionless_frequency -------------------------------------------
# This function calculates dimensionless frequencies, Q, Q', Qu, Qu', R
# Hsieh, 1987: Equations 13 & 14
# Rojstaczer and Riley, 1990
# Equations 4, 10, 19, 24, 28
.calc_dimensionless_frequency <- function(omega, thickness, diffusivity) {
  (thickness^2 * omega) / (2.0 * diffusivity)
}


# .calc_ohm ---------------------------------------------------------------
# Rojstaczer and Riley, 1990
# Equations 6, 23
# The specific storage term differs depending on the loading type
# Ss and Sa see Rojstaczer and Agnew, 1989 (Eq. 25) for discussion
# Sa is scaled by the term that accounts for one dimensional pore pressure
# diffusion on horizontal deformation which should be between 0.5 and 1.0.
# For stiff formations should be close to 1.0 Ss ~ Sa.  The difference
# therefore should be less than 2.0
.calc_ohm <- function(omega, specific_storage, k_vertical, specific_yield) {

  (1.0 - 1i) * sqrt((specific_storage * k_vertical) /
                      (2.0 * specific_yield^2 * omega))

}


# .calc_H1 ----------------------------------------------------------------
# Rojstaczer and Riley, 1990
# Equation 5
.calc_H1_H2 <- function(omega, ohm, thickness_aquifer, diffusivity_vertical) {


  # Rojstaczer and Riley, 1990 eq. 5ab
  # The equation in the paper appears to be missing a 2 in the calcuation for h1
  # as their equation does not reproduce the figures in the paper
  t1 <- 2.0 * (1i + 1.0) * sqrt(omega * thickness_aquifer^2 / (2 * diffusivity_vertical))

  # (minimize expensive operation)
  exp_t1   <- exp(t1)
  n_exp_t1 <- 1.0 / exp_t1

  # handle infinite values
  exp_t1   <- check_machine_max_complex(exp_t1)
  n_exp_t1 <- check_machine_max_complex(n_exp_t1)

  h1 <- 1.0 + n_exp_t1 - ohm * (1.0 - n_exp_t1)
  h2 <- 1.0 +   exp_t1 + ohm * (1.0 -   exp_t1)

  list(h1 = zapsmall(h1, 16), h2 = zapsmall(h2, 16))

}

check_machine_max_complex <- function(x, xmax = .Machine[['double.xmax']]) {
  complex(real      = ifelse(is.infinite(Re(x)), sign(Re(x)) * xmax, Re(x)),
          imaginary = ifelse(is.infinite(Im(x)), sign(Im(x)) * xmax, Im(x)))
}

check_machine_max <- function(x, xmax = .Machine[['double.xmax']]) {
  ifelse(is.infinite(x), sign(x) * xmax, x)
}


# .calc_u_v ---------------------------------------------------------------
# Rojstaczer and Riley, 1990
# Equation 27
.calc_u_v <- function(sqrt_Qu, h1, h2) {


  # (minimize expensive operation)
  cos_Qu <- cos(sqrt_Qu)
  sin_Qu <- sin(sqrt_Qu)

  # term 1 (minimize expensive operation)
  exp_t1   <- check_machine_max(exp(sqrt_Qu))
  n_exp_t1 <- zapsmall(1.0 / exp_t1, digits = 300)


  # Rojstaczer and Riley, 1990 eq. 27ab
  u <- n_exp_t1 * (-cos_Qu + sin_Qu) / (2.0 * sqrt_Qu * h1) +
    exp_t1 * ( cos_Qu + sin_Qu) / (2.0 * sqrt_Qu * h2) +
    ((1/h1) - zapsmall(1/h2, 16)) / (2.0 * sqrt_Qu)

  v <- n_exp_t1 * ( cos_Qu + sin_Qu) / (2.0 * sqrt_Qu * h1) +
    exp_t1 * (-cos_Qu + sin_Qu) / (2.0 * sqrt_Qu * h2) +
    (zapsmall(1.0/h2, 16) - (1/h1)) / (2.0 * sqrt_Qu)


  very_small <- 1e-300
  u <- ifelse(n_exp_t1 < very_small, 0.0, u)
  v <- ifelse(n_exp_t1 < very_small, 0.0, v)



  u <-  check_machine_max(u)
  v <-  check_machine_max(v)

  list(u = u, v = v)

}


# .calc_mn ----------------------------------------------------------------
# Rojstaczer, 1988 Determination of Fluid Flow Properties From the Response of
# Water Levels in Wells to Atmospheric Loading
# Equation 4
.calc_mn <- function(sqrt_r, attenuation) {

  m <- attenuation * (2 * cosh(sqrt_r) * cos(sqrt_r)) /
    (cosh(2 * sqrt_r) + cos(2 * sqrt_r))
  n <- attenuation * (2 * sinh(sqrt_r) * sin(sqrt_r)) /
    (cosh(2 * sqrt_r) + cos(2 * sqrt_r))

  # denom  <- (cosh(2.0 * (sqrt_r)) + cos(2.0 * (sqrt_r))) / (2.0 * attenuation)
  #
  # m      <- (cosh(sqrt_r) * cos(sqrt_r)) / denom
  # n      <- (sinh(sqrt_r) * sin(sqrt_r)) / denom

  # print(head(m))
  # print(head(n))

  m      <- ifelse(is.nan(m), 0.0, m)
  n      <- ifelse(is.nan(n), 0.0, n)

  list(m = m, n = n)

}



# .calc_h_w ---------------------------------------------------------------
# Cooper, Bredehoeft, Papadopulos, Bennett, 1965
# effective column height
# Equation 14, 15 (in text)
.calc_effective_height <- function(height_water, thickness_aquifer) {
  height_water + (3.0 * thickness_aquifer) / 8.0
}



#' convert_le_to_be
#'
#' @param tf complex transfer function
#'
#' @return the loading efficiency
#' @export
#'
convert_le_to_be <- function(tf) {

  (tf - 1)

}

#' convert_for_rojstaczer
#'
#' @param tf
#'
#' @return conjugate of tf
#' @export
#'
convert_for_rojstaczer <- function(tf) {

  Conj(tf)

}


# Hussein
# Determination of Fluid Flow Properties From the Response of Water Levels in Wells to Atmospheric Loading
#' areal_rojstaczer_semiconfined
#'
#'
#' @param frequency the frequency of the response
#' @param radius_well well radius
#' @param transmissivity aquifer transmissivity (L*L/t)
#' @param storage_confining confining layer storativity (L/L)
#' @param storage_aquifer aquifer storativity (L/L)
#' @param diffusivity_confining confining layer diffusivity
#' @param diffusivity_vadose air diffusivity of vadose zone
#' @param thickness_confining confining layer thickness
#' @param thickness_vadose vadose thickness
#' @param loading_efficiency the loading efficiency of the aquifer
#' @param attenuation an attenuation factor
#'
#' @return complex response vector in frequency domain
#' @export
#'
areal_rojstaczer_semiconfined <- function(frequency,
                                          radius_well,
                                          transmissivity,
                                          storage_confining,
                                          storage_aquifer,
                                          diffusivity_confining,
                                          diffusivity_vadose,
                                          thickness_confining,
                                          thickness_vadose,
                                          loading_efficiency,
                                          attenuation) {

  omega  <- .calc_omega(frequency)
  R      <- .calc_dimensionless_frequency(omega, thickness_vadose, diffusivity_vadose)
  Q      <- .calc_dimensionless_frequency(omega, thickness_confining, diffusivity_confining)
  W      <- .calc_W(omega, radius_well, transmissivity)

  sqrt_R <- sqrt(R)

  mn     <- .calc_mn(sqrt_R, attenuation)

  p0     <- ((mn$m - 1i * mn$n) - loading_efficiency) *
    exp(-(1i + 1.0) * sqrt(Q)) + loading_efficiency

  k0     <- Bessel::BesselK(
    ((W^2 * (storage_aquifer^2 + (storage_confining / (2.0 * Q))^2) )^0.25) *
      exp(0.5 * 1i * atan(2.0 * Q)), 0)  # well/aquifer

  #W <- 0#fifelse(Mod(W * k0) < 1e-10, 0.0, W)
  (-1 + p0) / (1.0 + (1i * 0.5 * W * k0))

  # return(data.table(frequency, R, Q, W, Q_div_W = Q / W, R_div_Q = R / Q, response = x0))
}


#' @inheritParams areal_rojstaczer_semiconfined
#' @export
areal_rojstaczer_unconfined <- function(frequency,
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
                                        attenuation) {


  omega <- .calc_omega(frequency)
  R     <- .calc_dimensionless_frequency(omega, thickness_vadose, diffusivity_vadose)
  Qu    <- .calc_dimensionless_frequency(omega, thickness_saturated_well, diffusivity_vertical)

  sqrt_R  <- sqrt(R)
  sqrt_Qu <- sqrt(Qu)

  mn      <- .calc_mn(sqrt_R, attenuation)
  ohm     <- .calc_ohm(omega, storage_aquifer, k_vertical, specific_yield)
  h1_h2   <- .calc_H1_H2(omega, 0, thickness_aquifer, diffusivity_vertical)
  uv      <- .calc_u_v(sqrt_Qu, h1_h2[['h1']], h1_h2[['h2']])

  # print(uv)
  p0 <- ((mn$m - 1i * mn$n) - loading_efficiency) *
    (uv$u + 1i * uv$v) + loading_efficiency


  (-1.0 + p0)


}
