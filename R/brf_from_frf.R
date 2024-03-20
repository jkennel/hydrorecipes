#' brf_from_frf
#'
#' Calculate the barometric response function from the transfer function
#'
#' @param x the complex transfer function vector
#'
#' @return vector of the barometric response function
#'
#' @export
#'
brf_from_frf <- function(x) {

  n <- floor(length(x) / 2) + 1

  imp <- fft(x, inverse = TRUE) / length(x)
  imp <- Mod(imp) * sign(Re(imp))

  cumsum(rev(imp)[1:n] + (imp)[1:n])

}
