test_that("kernel_apply works", {
  n <- 100L
  mdan <- modified_daniell(c(3, 3, 7), n)
  m <- matrix(rnorm(600), ncol = 6)
  a <- (kernel_apply(m, mdan))
  aa <- Re(a)[14:87]
  b <- kernel(coef = "modified.daniell", c(3, 3, 7))

  b <- kernapply(m[, 1], b)

  expect_equal(aa, b)

  expect_error(kernel_apply(m, modified_daniell(c(3, 3, 7), n + 1)))


})
