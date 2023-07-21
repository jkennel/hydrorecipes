#' be_storage
#'
#' This function estimates the storage coefficient from physical properties and the barometric efficiency
#' The equation used is from Batu 1998, eq 2-101, pg. 72
#'
#' @param be the barometric efficiency (numeric)
#' @param n the porosity 0 to 1 (numeric)
#' @param b aquifer thickness (m) (numeric) 
#' @param gamma the specific weight of water (N/m^3) (numeric)
#' @param beta the compressibility of water (4.786e-10 m^2/N) (numeric)
#'
#' @return storage calculated from barometric efficiency
#' 
#' @references Batu, V. (1998). Aquifer hydraulics: a comprehensive guide to 
#' hydrogeologic data analysis. John Wiley & Sons. 
#' 
#' @export
#'
#' @examples
#' be_storage(0.5, 0.32, 45)
#'
be_storage <- function(be,
                       n,
                       b,
                       gamma = 9799.74,
                       beta = 4.786e-10){
  
  return((gamma * n * beta * b) / be)
  
}
