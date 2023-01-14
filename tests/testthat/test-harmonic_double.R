test_that("harmonic_double works", {

  n <- 100L
  t <- seq(0, 86400 * 2, 3600)
  m <- matrix(1:(8*n), ncol = 8)
  cycles <- c(1,2,4,5)

  t_scale <- ((2*pi / 86400) * t-t[1]) %*% t(cycles)
  hd_r <- cbind(sin(t_scale), cos(t_scale))
  hd_c <- harmonic_double(t, cycles, 86400)
  expect_equal(hd_r, hd_c)

})

