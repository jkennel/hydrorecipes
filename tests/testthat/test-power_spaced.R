test_that("power_spaced works", {
  ps <- power_spaced(1000, 1, 10, 2)
  expect_equal(max(ps), 10)
  expect_equal(min(ps), 1)
  expect_equal(length(ps), 1000)
})
