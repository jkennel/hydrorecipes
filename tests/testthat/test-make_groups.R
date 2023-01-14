test_that("make_groups works", {
  n  <- 20034
  ng <- 200
  ma <- 2
  g  <- make_groups(n, power = 2, n_groups = ng, min_aggregate = ma)


  # even number
  expect_equal(g[1], 1)
  expect_equal(g[length(g)], 1)
  expect_equal(length(g), ng)
  expect_equal(min(g[-c(1, length(g))]), ma)
  expect_equal(sum(g), n)


  # odd number
  n <- 20031
  g <- make_groups(n, power = 2, n_groups = ng, min_aggregate = ma)
  expect_equal(g[1], 1)
  expect_equal(g[length(g)], 1)
  expect_equal(length(g), ng)
  expect_equal(min(g[-c(1, length(g))]), ma)
  expect_equal(sum(g), n)


  # same number groups as n
  n  <- 200
  ma <- 1
  ng <- 200
  expect_warning(g  <- make_groups(n, power = 2, n_groups = ng, min_aggregate = ma))
  expect_warning(g  <- make_groups(n, power = 2, n_groups = ng, min_aggregate = 3))
  expect_equal(max(g[1]), 1)
  expect_equal(min(g[1]), 1)
  expect_length(g, n)
  expect_equal(sum(g), n)


  # not enough values for the group specification
  n  <- 202
  ma <- 2
  expect_error(make_groups(n, power = 2, n_groups = ng, min_aggregate = ma))
  n  <- 400
  expect_error(make_groups(n, power = 2, n_groups = ng, min_aggregate = ma))

  # power and min_aggregate are too large
  n  <- 500
  ma <- 1
  expect_error(make_groups(n, power = 2, n_groups = ng, min_aggregate = ma))
  ma <- 3
  expect_error(make_groups(n, power = 2, n_groups = ng, min_aggregate = ma))

  ma <- 1
  expect_length(make_groups(n, power = 0.001, n_groups = ng, min_aggregate = ma), ng)
  expect_equal(sum(make_groups(n, power = 0.001, n_groups = ng, min_aggregate = ma)), n)

})
