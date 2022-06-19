test_that("lags_cpp works", {
  expect_error(check_lag(5, 10, 0))
  expect_error(distributed_lag_parallel(rnorm(10), matrix(1:20, ncol = 4), 5, 3, 4))
  expect_error(distributed_lag_parallel(rnorm(10), matrix(1:20, ncol = 4), 5, 3, -1))
  expect_error(distributed_lag_parallel(rnorm(10), matrix(1:20, ncol = 4), 7, 5, 4))

})
