test_that("multiply_ffts works", {

  m <- matrix(rnorm(9e2), ncol = 3)

  a <- fft_matrix(m, n_new = nrow(m))
  b <- mvfft(m)

  d <- multiply_ffts(a)

  e <- b * Conj(b)

  # check non-cross terms
  expect_equal(d[,1], e[,1])
  expect_equal(d[,4], e[,2])
  expect_equal(d[,6], e[,3])


  })
