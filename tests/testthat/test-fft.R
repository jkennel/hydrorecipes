test_that("fft_matrix works", {

  m <- matrix(rnorm(300), ncol = 3)
  a <- fft_matrix(m, n_new = nrow(m))
  b <- mvfft(m)
  expect_equal(a, b)

})

test_that("convolve_vec works", {

  a <- convolve_vec(1:100, as.matrix(1:100))

  b <- rev(convolve(1:100, rev(1:100), type = 'circular'))

  expect_equal(a, b)

})

test_that("convolve_matrix works", {

  a <- convolve_matrix(x = 1:100,
                       y = as.matrix(1:10),
                       remove_partial = TRUE,
                       reverse = TRUE)

  b <- convolve(1:100, rev(1:10), type = 'filter')


  expect_equal(tail(b, 91), tail(a[, 1], 91))

})


test_that("spec_pgram works", {
  tap <- 0.1
  u2 <- (1.0 - (5/8) * tap * 2.0)

  n  <- 2000
  nc <- 2
  n_out <- n/nc/2
  m <- matrix(rnorm(n), ncol = nc)

  a <- Re(spec_pgram(m,
                     spans = 3,
                     detrend = FALSE,
                     demean = FALSE,
                     taper = tap))

  a <- a[-1,] / u2
  b <- spec.pgram(m,
                  detrend = FALSE,
                  demean = FALSE,
                  plot = FALSE,
                  spans = 3,
                  taper = tap)$spec

  expect_equal(a[1:n_out, c(1,3)], b)


  m2 <- sweep(m, 2, colMeans(m), check.margin=FALSE)
  a <- Re(spec_pgram(m2,
                     spans = 3,
                     detrend = FALSE,
                     demean = FALSE,
                     taper = tap))
  a <- a[-1,] / u2
  b <- spec.pgram(m,
                  detrend = FALSE,
                  demean = TRUE,
                  plot = FALSE,
                  spans = 3,
                  taper = tap)$spec

  expect_equal(a[1:n_out, c(1,3)], b)


  a <- Re(spec_pgram(m,
                     spans = 3,
                     detrend = FALSE,
                     demean = TRUE,
                     taper = tap))
  a <- a[-1,] / u2
  b <- spec.pgram(m, detrend = FALSE,
                  demean = TRUE, plot = FALSE, spans = 3, taper = tap)$spec

  expect_equal(a[1:n_out, c(1,3)], b)


  a <- Re(spec_pgram(m,
                     spans = 3,
                     detrend = TRUE,
                     demean = FALSE,
                     taper = tap))
  a <- a[-1,] / u2
  b <- spec.pgram(m, detrend = TRUE, demean = FALSE, plot = FALSE, spans = 3,
                  taper = tap)$spec

  expect_equal(a[1:n_out, c(1,3)], b)

})

test_that("spec_welch works", {

  scale <- 0.37
  x <- rnorm(10000)
  y <- x * scale
  m <- matrix(c(y,x), ncol = 2)

  a <- spec_welch(m,
                  length_subset = nrow(m)/1.005,
                  overlap = 0.995,
                  window = window_hann(nrow(m)/1.005))


  cp <- ordinary_coherence_phase(a)
  expect_equal(min(cp[,1]), 1)
  expect_equal(max(cp[,1]), 1)
  expect_equal(min(cp[,2]), 0)
  expect_equal(max(cp[,2]), 0)

  solve_cplx_irr(a,
                 power = 2,
                 n_groups = 100,
                 min_aggregate = 2)

})

test_that("ordinary_coherence_phase works", {

  scale <- 0.37
  x <- rnorm(1000)
  y <- x * scale
  m <- matrix(c(y,x), ncol = 2)
  tap <- 0.01

  a <- spec_pgram(m,
                  spans = 3,
                  detrend = FALSE,
                  demean = FALSE,
                  taper = tap)
  # coherence of 1
  # phase of 0
  cp <- ordinary_coherence_phase(a)
  expect_equal(min(cp[,1]), 1)
  expect_equal(max(cp[,1]), 1)

  expect_equal(min(cp[,2]), 0)
  expect_equal(max(cp[,2]), 0)


  set.seed(123)
  x <- rnorm(1000)
  m <- matrix(x, ncol = 2)

  a <- spec_pgram(m,
                  spans = c(3, 5, 9),
                  detrend = FALSE,
                  demean = FALSE,
                  taper = tap)

  cp <- ordinary_coherence_phase(a)
  # plot(cp[,1],type= 'l')
  expect_lt(max(abs(cp[,1])), 0.5)


})

test_that("transfer_function works", {

  scale <- 0.37
  x <- rnorm(1000)
  y <- x * scale
  m <- matrix(c(y, x), ncol = 2)
  tap <- 0.01

  pg <- spec_pgram(m, spans = 3, TRUE, TRUE, taper = tap)
  tf <- solve_cplx_irr(pg, 2, 100, min_aggregate = 1)
  expect_equal(range(Mod(tf)) - scale, rep(0.0, 2))
  expect_equal(range(Arg(tf)), rep(0.0, 2))


  pg <- spec_welch(m, length_subset = nrow(m)/1.005,
                   overlap = 0.995,
                   window = window_hann(nrow(m)/1.005))
  tf <- solve_cplx_irr(pg, 2, 100, min_aggregate = 1)
  expect_equal(range(Mod(tf)) - scale, rep(0.0, 2))
  expect_equal(range(Arg(tf)), rep(0.0, 2))


  scale <- 0.37
  x <- rnorm(1e4)
  y <- x * scale
  m <- matrix(c(y,x), ncol = 2)
  tf <- transfer_pgram(m,
                       spans = 3,
                       detrend = FALSE,
                       demean = FALSE,
                       0.1)

  expect_equal(range(Mod(tf)) - scale, rep(0.0, 2))
  expect_equal(range(Arg(tf)), rep(0.0, 2))
  tf <- transfer_welch(m,
                       length_subset = as.integer(nrow(m) / 1.1),
                       overlap = 0.9,
                       window = window_hann(as.integer(nrow(m) / 1.1)))
  expect_equal(range(Mod(tf)) - scale, rep(0.0, 2))
  expect_equal(range(Arg(tf)), rep(0.0, 2))
  # plot(Mod(tf))
})


test_that("determine_frequency works", {
  n <- 1000
  expect_equal(range(determine_frequency(n)), c(0, 0.5-1.0/n))
})


test_that("group_frequency works", {

  n  <- 1000
  g  <- determine_frequency(n)
  gf <- group_frequency(g, 2, 30, 1)

  expect_equal(range(gf), range(g))
  expect_equal(length(gf), 30)

})

# test_that("frf_to_brf works", {
#   library(aquifer)
#   library(data.table)
#   data("transducer")
#   setDT(transducer)
#
#   max_syn_lag <- 86400/120
#   # vadose kernel -----------------------------------------------------------
#   kern_vad_a <- aquifer::vadose_response(time = as.numeric(0:max_syn_lag),
#                                          D = 0.1,
#                                          L = 40,
#                                          precision = 1e-16,
#                                          inverse = TRUE)
#
#   kern_vad_b <- aquifer::vadose_response(time = as.numeric(0:max_syn_lag),
#                                          D = 0.5,
#                                          L = 40,
#                                          precision = 1e-16,
#                                          inverse = TRUE)
#   kern_vad_comb <- (kern_vad_a+kern_vad_b) / 2
#
#   # remove small values
#   kern_vad_comb[kern_vad_comb < 8e-16] <- 0
#   kern_vad_a[kern_vad_a < 8e-16] <- 0
#
#
#
#   # slug kernel -------------------------------------------------------------
#   t <- c(1e-6, 1:max_syn_lag)
#   kern_slug <- c(slug_cbp(t,
#                           radius_screen = 0.10,
#                           storativity = 1e-5,
#                           transmissivity = 5e-4,
#                           n = 8L))
#   s <- max(kern_slug)
#
#   # combined kernel ---------------------------------------------------------
#   kern_comb <- as.matrix(rev(c(0, diff(kern_vad_a))) + (c(s, diff(kern_slug))))
#   plot(1- (cumsum(kern_comb)), type = 'l', log = 'x')
#   x <- as.matrix(cumsum(rnorm(86400 * 3)))
#
#   # transducer[, wl_comb := waterlevel::synthetic_wl(baro, datetime, kernel = rev(kern_comb))]
#
#
#   y <- convolve_matrix(x,
#                        kern_comb[length(kern_comb):1,, drop = FALSE],
#                        remove_partial = FALSE,
#                        reverse = TRUE)
#   plot(x, type = 'l')
#   points(-y, type = 'l', col = 'red')
#
#   m <- cbind(x, y)
#   tap <- 0.01
#   pg <- spec_pgram_trunc(m, spans = c(3), TRUE, TRUE, taper = tap)
#   tf <- solve_cplx_irr(pg, 2, 50, min_aggregate = 1)
#
#   plot(Mod(tf[,1]), type = 'h', log = 'xy')
# })



# test_that("interpolate_tf works", {
#   scale <- 0.37
#   x <- rnorm(1e4)
#   y <- x * scale
#   m <- matrix(c(y,x), ncol = 2)
#   tf <- transfer_pgram(m, spans = 3, detrend = FALSE, demean = FALSE,
#                              0.1, 2, 50)
#   g  <- make_groups(1e4, power = 2, n_groups = 50, min_aggregate = 2)
#   at <- cumsum(g)-1
#   # aa <- interpolate_tf(at,
#   #                       y = tf,
#   #                       knots = c(0, 5, 10, 50, 9999),
#   #                       degree = 3)
#
#   expect_equal(range(cumsum(aa) - scale), rep(0.0, 2))
#
#
# })
