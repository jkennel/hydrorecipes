test_that("bs_eigen equal knot range works", {

  # equal knot range to inputs
  hy_bs <- b_spline(0:100, knots = c(0, 10, 50, 100), degree = 3)
  sp_bs <- splines::bs(0:100, knots = c(10, 50), degree = 3, intercept = TRUE)
  expect_equal(hy_bs, sp_bs[, 1:6], ignore_attr = TRUE)
})


test_that("bs_eigen no intercept works", {

  # no intercept
  hy_bs <- b_spline(0:100, knots = c(0, 10, 50, 100), degree = 3)[,-1]
  sp_bs <- splines::bs(0:100, knots = c(10, 50), degree = 3, intercept = FALSE)
  expect_equal(hy_bs, sp_bs[, 1:5], ignore_attr = TRUE)
})


test_that("bs_eigen larger knot range works", {

  # larger knot range than inputs
  hy_bs <- b_spline(0:100, knots = c(-1, 10, 50, 101), degree = 3)
  sp_bs <- splines::bs(0:100, knots = c(10, 50), Boundary.knots = c(-1, 101), degree = 3, intercept = TRUE)
  expect_equal(hy_bs, sp_bs[, 1:6], ignore_attr = TRUE)
})


test_that("bs_eigen smaller knot range works", {

  # check which indices (extrapolate)
  hy_bs <- b_spline(0:100, knots = c(2, 10, 50, 99), degree = 3)
  suppressWarnings(sp_bs <- splines::bs(0:100,
                                        knots = c(10, 50),
                                        Boundary.knots = c(2, 99),
                                        degree = 3, intercept = TRUE))
  expect_equal(hy_bs, sp_bs[, 1:6], ignore_attr = TRUE)

})

test_that("bs_eigen different degree works", {
  # degree == 4
  hy_bs <- b_spline(0:100, knots = c(0, 10, 50, 100), degree = 4)
  sp_bs <- splines::bs(0:100, knots = c(10, 50),
                       Boundary.knots = c(0, 100),
                       degree = 4, intercept = TRUE)
  expect_equal(hy_bs, sp_bs[, 1:7], ignore_attr = TRUE)

  # degree == q
  hy_bs <- b_spline(0:1e2, knots = c(0, 10, 50, 100), degree = 2)
  sp_bs <- splines::bs(0:1e2, knots = c(10, 50),
                       Boundary.knots = c(0, 100),
                       degree = 2, intercept = TRUE)

  expect_equal(hy_bs, sp_bs[, 1:5], ignore_attr = TRUE)
})
