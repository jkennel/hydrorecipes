test_that("log_lags works", {

  expect_equal(min(log_lags(10, 100)), 0)
  expect_equal(max(log_lags(10, 100)), 100)
  expect_equal(log_lags(1, 100), 0)
  expect_equal(max(log_lags(2, 100)), 100)
  expect_equal(length(log_lags(1, 100)), 1)
  expect_equal(length(log_lags(101, 100)), 101)
  expect_equal(log_lags(1, 0), 0)

  expect_warning(log_lags(200, 100))

  expect_error(log_lags(1, -1))
  expect_error(log_lags(0, 100))
  expect_error(log_lags(-1, 100))
  expect_error(log_lags('a', 100))
  expect_error(log_lags(1, 's'))
  expect_error(log_lags('a', 's'))
  expect_error(log_lags(c(1,2), 100))
  expect_error(log_lags(1, c(20, 21)))

})
