test_that("spec_welch padding works", {
  n_var <- 4
  x <- matrix(rnorm(1e4 * n_var), ncol = n_var)

  # check number of rows of output
  out <- spec_welch(x, length_subset = 255, overlap = 0.5, window = window_hann(255))
  expect_equal(nrow(out), nextn(256))
  expect_equal(ncol(out), n_var * (n_var +1)/2)



  # m <- rnorm(10000)
  # mm <- as.matrix(m)
  # microbenchmark::microbenchmark(
  #   a <- oce::pwelch(m, fs = 100, plot = FALSE),
  #   b <- spec_pgram(mm, spans = 3,taper = 0.1, detrend = TRUE, demean = FALSE)
  # )

})

test_that("spec_welch one column matrix works", {
  x <- matrix(rnorm(1e3), ncol = 1)
  expect_silent(spec_welch(x, length_subset = 255, overlap = 0, window = window_hann(255)))
})


test_that("spec_welch input checks works", {
  x <- matrix(rnorm(1e3 * 4), ncol = 4)

  # check inputs
  expect_error(spec_welch_trunc(x, length_subset = 55, overlap = 1.1, window = window_hann(255)))
  expect_error(spec_welch_trunc(x, length_subset = 55, overlap = -1, window = window_hann(255)))
  expect_error(spec_welch_trunc(x, length_subset = 1001, overlap = 0.5, window = window_hann(255)))

})


# test_that("spec_welch spec works", {
#   library(oce)
#   tap <- 0.1
#   u2 <- (1.0 - (5/8) * tap * 2.0)
#
#   Fs <- 1e6
#   t <- seq(0, 0.296, 1/Fs)
#   x <- cos(2 * pi * t * 200) + rnorm(n=length(t))
#   X <- ts(x, frequency=Fs)
#
#   nfft <- 2^8
#
#   a <- Re(spec_welch(as.matrix(X),
#                               length_subset = nfft,
#                               overlap = 0,
#                               window = hydrorecipes:::window_rectangle(nfft)))
#   b <- oce::pwelch(X,
#                    nfft = nfft,
#                    noverlap = 0,
#                    detrend = FALSE,
#                    demean = FALSE,
#                    plot = FALSE,
#                    window = hydrorecipes:::window_rectangle(nfft))$spec
#
#   microbenchmark::microbenchmark(
#     a <- Re(spec_welch(as.matrix(X),
#                                 length_subset = nfft,
#                                 overlap = 0,
#                                 window = hydrorecipes:::window_rectangle(nfft))),
#     b <- oce::pwelch(X,
#                      nfft = nfft,
#                      noverlap = 0,
#                      detrend = FALSE,
#                      demean = FALSE,
#                      plot = FALSE,
#                      window = hydrorecipes:::window_rectangle(nfft))$spec,
#     d <- ravetools::pwelch(X,fs = 1,
#                            window = nfft,
#                            nfft = 1,
#                      noverlap = 0),
#     times = 1
#   )
#   a <- a[-1,] / Fs
#   d <- d$spec
#   a/d
#
#   plot(b, type = 'l')
#   points(d[1:64], type = 'l')
#   plot(d[1:64], type = 'l')
#   points(a[1:(nfft/2)], type = 'l', col = 'red')
#   expect_equal(a[1:64], b)
# })
