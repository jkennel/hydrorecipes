test_that("lags_cpp works", {
  expect_error(check_lag(5, 10, 0))
  expect_error(distributed_lag_parallel(rnorm(10), matrix(1:20, ncol = 5), 5, 3, 4))
  expect_error(distributed_lag_parallel(rnorm(10), matrix(1:20, ncol = 5), 5, 3, -1))
  expect_error(distributed_lag_parallel(rnorm(10), matrix(1:20, ncol = 5), 7, 5, 4))
  expect_error(distributed_lag_parallel(rnorm(10), matrix(1:20, ncol = 5), 4, 0, 0))


  expect_equal(nrow(distributed_lag_parallel(rnorm(10), matrix(1:20, ncol = 5), 4, 2, 0)), 5)
  expect_equal(nrow(distributed_lag_parallel(rnorm(10), matrix(1:20, ncol = 5), 4, 3, 0)), 4)


  expect_error(shift_subset(rnorm(10), n_subset = 30, lag = 10))
  expect_error(shift_subset(rnorm(10), n_subset = 30, lag = 1))
  expect_error(shift_subset(rnorm(10), n_shift = 30, lag = 1))



})
