test_that("determine_frequency works", {
  n <- 100
  f <- determine_frequency(n)
  expect_equal(max(f), 0.5 - 1/n)
  expect_equal(min(f), 0)
  expect_equal(length(f), n/2)
})
