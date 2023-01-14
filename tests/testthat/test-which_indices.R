test_that("which_indices works", {

  # first and last knots are the end points.  The endpoints are extended if the
  # knots aren't the full range
  wh_1 <- which_indices(0:100, c(2, 10, 50, 99))
  wh_2 <- which_indices(0:100, c(0, 10, 50, 100))

  expect_equal(wh_1, wh_2)

})
