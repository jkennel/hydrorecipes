#' barker_herbert_impulse_2 <- function(p,
#'                                      flow_rate,
#'                                      radius,
#'                                      radius_patch,
#'                                      t_1,
#'                                      t_2,
#'                                      s_1,
#'                                      s_2) {
#'
#'   N  <- sqrt(s_1 * p / t_1)
#'   A  <- sqrt(s_2 * p / t_2)
#'   ct <- (t_2 / t_1) * (A / N)
#'   s  <- flow_rate / (2.0 * pi * t_1)
#'
#'   bi_n <- Bessel::BesselI(N * radius_patch, 0L, nSeq = 2L)
#'   bk_a <- Bessel::BesselK(A * radius_patch, 0L, nSeq = 2L)
#'   bk_n <- Bessel::BesselK(N * radius_patch, 0L, nSeq = 2L)
#'
#'   denom <- (ct * bi_n[, 1L] * bk_a[, 2L] +
#'               bi_n[, 2L] * bk_a[, 1L]) * p
#'
#'   term_1 <- (bk_n[, 2L] * bk_a[, 1L]  -
#'                bk_a[, 2L] * bk_n[, 1L] * ct)
#'   term_2 <- (bk_n[, 1L] * bi_n[, 2L] +
#'                bk_n[, 2L] * bi_n[, 1L])
#'   # print(term_1)
#'   # print(term_2)
#'   drawdown <- matrix(NA, ncol = ncol(p) * length(radius), nrow = nrow(p))
#'
#'   for(i in seq_along(radius)) {
#'     start = 1L + ncol(p) * (i-1L)
#'     end = ncol(p) * i
#'     d <- radius[i]
#'     if (d <= radius_patch) {
#'       # print(Bessel::BesselK(N * d, 0) / p)
#'       # print(((term_1 * Bessel::BesselI(N * d, 0)) / denom))
#'       drawdown[, start:end] <- (Bessel::BesselK(N * d, 0) / p  +
#'                                   (term_1 * Bessel::BesselI(N * d, 0)) / denom)
#'
#'     } else {
#'       drawdown[, start:end] <- (term_2 * Bessel::BesselK(A * d, 0) / denom)
#'     }
#'   }
#'
#'   (s * drawdown)
#' }
#'
#'
#' #' barker_herbert
#' #'
#' #' @param times
#' #' @param distance
#' #' @param radius
#' #' @param q
#' #' @param t_1
#' #' @param t_2
#' #' @param s_1
#' #' @param s_2
#' #'
#' #' @return drawdown using barker and herbert 1982
#' #' @export
#' #'
#' barker_herbert <- function(time,
#'                            radius,
#'                            radius_patch,
#'                            flow_rate,
#'                            t_1,
#'                            t_2,
#'                            s_1,
#'                            s_2) {
#'
#'   stehfest(time, n = 12L, impulse = barker_herbert_impulse_2, flow_rate, radius, radius_patch,
#'            t_1, t_2, s_1, s_2)
#'
#' }
#'
#'
