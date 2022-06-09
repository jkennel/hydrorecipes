#https://en.wikipedia.org/wiki/Window_function

# ==============================================================================
#' @title
#' fftw_convolve
#'
#' @description
#' Do convolution for a data series and matrix of equal length filters.
#' For long data series this function will be much faster than non-fft methods.
#'
#' @param x \code{numeric} vector or single column matrix the data series
#' @param y \code{numeric} matrix of filter(s)
#' @param normalize \code{logical} do you want the values normalized
#' @param align \code{character} what alignment for the convolution
#'
#' @return convolution of data series and filter
#'
#' @export
#'
#' @examples
#'
#' library(splines)
#'
#' knots <- round(10^seq(0.5, 4.635484, length.out = 9))
#' n_lags <- 86400*2
#'
#' bla <- ns(0:(n_lags), knots = knots,
#'           Boundary.knots = c(0, n_lags))
#' bv <- rnorm(86400*2 + 3600*24*4)
#' system.time({
#' dl <- fftw_convolve(bv, bla, normalize = FALSE)
#' })
#'
fftw_convolve <- function(x,
                          y,
                          normalize = TRUE,
                          align = 'center') {

  y <- as.matrix(y)
  x <- as.matrix(x)

  n_x_in <- nrow(x)
  n_y    <- nrow(y)

  # reverse order
  y <- y[n_y:1, , drop = FALSE]

  if (n_y > n_x_in) {
    stop('The length of y should be less than or equal to x')
  }


  # set NA indices
  start <- 1:(ceiling(n_y / 2) - 1)
  end   <- (n_x_in - floor(n_y / 2) + 1):(n_x_in)


  if(n_x_in %% 2 == 1) {
    x     <- c(x, 0.0)
    n_x   <- length(x)
    x_pad <- n_x / 2
    sub   <- c((n_x + x_pad + 1):(n_x + 2 * x_pad), 1:(x_pad-1))
  } else {
    n_x   <- n_x_in
    x_pad <- n_x / 2
    sub   <- c((n_x + x_pad + 1):(n_x + 2 * x_pad), 1:(x_pad))
  }

  n_sub <- length(sub)

  # normalize filter to sum to 1
  if (normalize) {
    y <- t(t(y)  / colSums(y))
  }

  y_pad <- x_pad + ceiling((n_x - n_y) / 2)

  x_pad <- rep(0.0, x_pad)
  y_pad <- rep(0.0, y_pad)

  # do the FFT on the x vector
  f <- fftw::FFT(c(x_pad, x, x_pad))
  u <- matrix(NA_real_, nrow = n_sub, ncol = ncol(y))

  if (n_y %% 2 == 0) {

    warning('Values are shifted 0.5 units forward. Use odd number filter for better centering')
    for (i in 1:ncol(y)){
      u[, i] <- Re(fftw::IFFT(fftw::FFT(c(y_pad, y[,i], y_pad)) * f))[sub]
    }

  } else {
    for (i in 1:ncol(y)){
      u[,i] <- Re(fftw::IFFT(fftw::FFT(c(y_pad, y[,i], y_pad[-1])) * f))[sub]
    }

  }

  u[start,] <- NA_real_
  u[end,]   <- NA_real_

  # adjust output for alignment
  if (align == 'right') {
    for(i in 1:ncol(y)){
      u[,i] <- data.table::shift(u[,i], n = length(end), type = 'lag')
    }
  } else if (align == 'left') {
    for(i in 1:ncol(y)) {
      u[,i] <- data.table::shift(u[,i], n = length(start), type = 'lead')
    }
  }

  return(u)
}
#' window_rectangular
#'
#' Rectangular window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_rectangular(100)
#'
window_rectangular <- function(n){

  rep(1.0, n)

}


#' window_hann
#'
#' Hann window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_hann(100)
#'
window_hann <- function(n){

  0.5 * (1 - cos((2 * pi * 0:(n - 1)) / (n - 1)))

}



#' window_hamming
#'
#' Hamming window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_hamming(100)
#'
window_hamming <- function(n){

  0.54 - (0.46 * cos((2 * pi * 0:(n - 1)) / (n - 1)))

}

#' window_blackman
#'
#' Blackman window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_blackman(100)
#'
window_blackman <- function(n){

  a0 <- 0.42
  a1 <- 0.5
  a2 <- 0.08
  k <- 0:(n - 1) / (n - 1)

  a0 - a1 * cos(2 * pi * k) +
    a2 * cos(4 * pi * k)

}

#' window_first_deriv
#'
#' First derivative window for FFT
#'
#' @param n length of signal (integer)
#' @param a0 \code{numeric} coefficient
#' @param a1 \code{numeric} coefficient
#' @param a2 \code{numeric} coefficient
#' @param a3 \code{numeric} coefficient
#'
#' @return window
#'
#' @export
#'
#' @examples
#' # nuttall window
#' window_first_deriv(100, 0.355768, 0.487396, 0.144232, 0.012604)
#'
window_first_deriv <- function(n,
                               a0,
                               a1,
                               a2,
                               a3) {

  k <- 0:(n - 1) / (n - 1)
  a0 - a1 * cos(2 * pi * k) +
    a2 * cos(4 * pi * k) -
    a3 * cos(6 * pi * k)

}

#' window_nuttall
#'
#' Nuttall window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_nuttall(100)
#'
window_nuttall <- function(n){

  a0 <- 0.355768
  a1 <- 0.487396
  a2 <- 0.144232
  a3 <- 0.012604

  window_first_deriv(n, a0, a1, a2, a3)

}

#' window_blackman_nuttall
#'
#' Blackman-Nuttall window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_blackman_nuttall(100)
#'
window_blackman_nuttall <- function(n){

  a0 <- 0.3635819
  a1 <- 0.4891775
  a2 <- 0.1365995
  a3 <- 0.0106411

  window_first_deriv(n, a0, a1, a2, a3)

}

#' window_blackman_harris
#'
#' Blackman-Harris window for FFT
#'
#' @param n length of signal (integer)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_blackman_harris(100)
#'
window_blackman_harris <- function(n){

  a0 <- 0.35875
  a1 <- 0.48829
  a2 <- 0.14128
  a3 <- 0.01168

  window_first_deriv(n, a0, a1, a2, a3)

}



#' window_gaussian
#'
#' Gaussian window for FFT
#' https://en.wikipedia.org/wiki/Window_function
#'
#' @param n length of signal (integer)
#' @param sigma the standard deviation in periods (numeric)
#'
#' @return window
#'
#' @export
#'
#' @examples
#' window_gaussian(100, 0.3)
#
window_gaussian <- function(n,
                            sigma) {

  exp(-0.5 * ((0:(n-1) - (n-1) / 2) / (sigma * (n-1) / 2))^2)

}



