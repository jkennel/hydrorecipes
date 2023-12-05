#' #' stehfest
#' #'
#' #' @description Inverse Laplace transform using Stehfest algorithm
#' #'
#' #' @param time elapsed time since slug introduction
#' #' @param n coefficient for Stehfest algorithm (should be < 20)
#' #' @param impulse function to invert
#' #' @param ... arguments to pass to impulse
#' #'
#' #' @return inverse laplace transform of impulse function
#' #' @export
#' #'
#' stehfest <- function(time, n = 12L, impulse, ...) {
#'
#'   n_div_2 <- n / 2L
#'   ind     <- 1L:n
#'   fac     <- factorial(0L:n)
#'   v       <- rep(NA_real_, n)
#'
#'   for (i in ind) {
#'
#'     k <- min(i, n_div_2):trunc((i + 1L) / 2L)
#'     # print('-----------------')
#'     # print(n_div_2)
#'     # print(2L * k + 1L)
#'     # print(n_div_2 - k + 1L)
#'     # print(k + 1L)
#'     # print(k)
#'     # print(i - k + 1L)
#'     # print(2L * k - i + 1L)
#'
#'     z <- sum(((k^n_div_2) * fac[2L * k + 1L]) /
#'                (fac[n_div_2 - k + 1L] * fac[k + 1L] * fac[k] *
#'                   fac[i - k + 1L] * fac[2L * k - i + 1L]))
#'
#'     v[i] <- (-1)^(n_div_2 + i) * z;
#'
#'   }
#'   # print(v)
#'   time <- log(2.0) / time
#'   # print("time")
#'   # print(time)
#'   p    <- tcrossprod(as.numeric(ind), time)
#'   # print("p")
#'   # print(p)
#'   # print("impulse")
#'   # print(as.vector(crossprod(impulse(p, ...), v)))
#'   out  <- Re(as.vector(crossprod(impulse(p, ...), v) * time))
#'
#'   return(out)
#' }
