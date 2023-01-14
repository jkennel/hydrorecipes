test_that("pad_vector works", {

  x <- rnorm(10)

  expect_equal(pad_vector(x, 10, 20), c(x, rep(0, 10)))
  expect_equal(pad_vector(x, 10, 10), x)

  expect_error(pad_vector(x, 10, 9))

})
