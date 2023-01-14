test_that("window_hann works", {

  window_hann_r <- function(n){ 0.5 * (1 - cos((2 * pi * 0:(n - 1)) / (n - 1)))}

  expect_equal(window_hann(100), window_hann_r(100))
  expect_equal(window_hann(1), 1)
  expect_equal(window_hann(2), c(0.5, 0.5))

  expect_error(window_hann(0))


  expect_equal(window_hann(100), Re(window_hann_cplx(100)))
  expect_equal(window_hann(100), Im(window_hann_cplx(100)))



})

test_that("window_tukey works", {


  expect_equal(window_tukey(100, 1), window_hann(100))
  expect_equal(window_tukey(100, 0), rep(1.0, 100))


  expect_equal(window_tukey(100, 0.5), gsignal::tukeywin(100, 0.5))
  expect_equal(window_tukey(100, 0.25), gsignal::tukeywin(100, 0.25))
  expect_equal(window_tukey(101, 0.25), gsignal::tukeywin(101, 0.25))


  expect_equal(window_tukey(1, 0.25), gsignal::tukeywin(1, 0.25))
  expect_equal(window_tukey(2, 0.25), gsignal::tukeywin(2, 0.25))
  expect_equal(window_tukey(101, 0.25), gsignal::tukeywin(101, 0.25))


  expect_error(window_tukey(0))
})

# bspec gives different values than gsignal
# test_that("window_tukey works", {
#   a <- window_tukey(100, 0.2)
#   b <- bspec::tukeywindow(100, 0.2)
#   expect_equal(a, b)
#
#   a <- window_tukey(101, 0.1)
#   b <- bspec::tukeywindow(101, 0.1)
#   expect_equal(a, b)
# })

# need better tests here
test_that("window_scale works", {
  expect_gt(window_scale(window_hann(100), 100, 100),
            window_scale(window_hann(100), 100, 150))

})
