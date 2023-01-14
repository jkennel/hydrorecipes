#'
#' brf_from_frf <- function(x, dc1, dc2) {
#'
#'   # half-spectrum length
#'   n  <- length(x)
#'
#'   # reconstruct the full spectrum from half spectra
#'   x <- c(dc1, x, dc2, Conj(rev(x)))
#'   imp <- fftw::IFFT(x, scale = TRUE)
#'
#'   Re(imp[1:n] + imp[length(x):(length(x) - n + 1)])
#'
#' }
#'
#' # x are the indices
#' # y is a matrix of the complex frequency response
#' #' frequency_to_time_domain
#' #'
#' #' @param x values for the center of the complex response (numeric vector)
#' #' @param y values for the frequency response (complex vector)
#' #'
#' #' @return
#' #' @export
#' #'
#' #' @examples
#' frequency_to_time_domain <- function(x, y) {
#'
#'   dc1 <- y[1,] # DC values
#'   dc2 <- y[nrow(y),] # DC values
#'   x_n <- x[-c(1, length(x))] # DC indices
#'   y_n <- y[-c(1, nrow(y)),] # DC values
#'   knots <- log_lags(20, max(x_n))
#'   knots <- knots[-c(1, length(knots))]
#'   len <- min(knots):max(knots)
#'
#'   # these are used to smooth and interpolate the complex response
#'   sp_in  <- splines2::bSpline(x_n, knots = knots, Boundary.knots = c(0, max(x_n)+1e-16))
#'   sp_out <- splines2::bSpline(len, knots = knots, Boundary.knots = c(0, max(x_n)+1e-16))
#'
#'
#'   out <- matrix(NA_real_,
#'                 nrow = length(len),
#'                 ncol = ncol(y_n) + 1)
#'   out[, 1] <- len
#'
#'   # loop through each transfer function - smooth the response - estimate brf
#'   for (i in 1:ncol(y_n)) {
#'
#'     # smooth the real and imaginary components separately
#'     re <- Re(y_n[, i])
#'     im <- Im(y_n[, i])
#'
#'     fit_r <- RcppEigen::fastLm(X = sp_in, y = re)$coefficients
#'     fit_i <- RcppEigen::fastLm(X = sp_in, y = im)$coefficients
#'
#'     re <- as.vector(sp_out %*% fit_r)
#'     im <- as.vector(sp_out %*% fit_i)
#'
#'     # estimate the brf from the smoothed uniform spaced frf
#'     # out[, i + 1] <- hydrorecipes:::frf_to_brf(complex(real = re, imaginary = im),
#'     #                                           dc1[1], dc2[1])
#'
#'     out[, i + 1] <- hydrorecipes:::brf_from_frf(complex(real = re, imaginary = im),
#'                                                 dc1[1], dc2[1])
#'
#'   }
#'
#'
#'   out
#' }
