#' baro_elev
#'
#' Calculate the effect of elevation on barometric pressure.
#'
#' @param elevation  vector elevations (numeric)
#' @param pressure_sea standard pressure (numeric)
#' @param temperature  temperature in celcius (numeric)
#' @param gravity  gravity (numeric)
#'
#' @return barometric pressure dependence on elevation
#'
#' @export
#'
#' @examples
#' plot(baro_elev(0:10000), type='l', xlab = "elev (m)", ylab = 'pressure')
#'
baro_elev <- function(elevation = 0.0,        # meters (m)
                      pressure_sea = 101325.0,  # Pascal (Pa)
                      temperature = 20.0,     # celcius (C)
                      gravity = 9.80665) {    
  
  gas_const <- 8.31432
  molar_mass_air <- 0.0289644
  
  pressure_sea * exp(-gravity * molar_mass_air * (elevation) / (gas_const * (temperature + 273.15)))
  
}