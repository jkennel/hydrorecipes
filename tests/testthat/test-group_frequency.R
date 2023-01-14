test_that("group_frequency works", {

  n <- 2034
  f <- 1:n
  ng <- 100
  ma <- 2
  g  <- group_frequency(f, power = 3, n_groups = ng, min_aggregate = ma)

  expect_equal(length(g), ng)
  expect_equal(min(g), min(f))
  expect_equal(max(g), max(f))

})

