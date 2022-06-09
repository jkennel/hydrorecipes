
# log_lags ----------------------------------------------------------------


expect_silent(log_lags(1L, 20L))

expect_silent(log_lags(5, 5))

expect_warning(log_lags(10, 5))

expect_error(log_lags(0, 20))

expect_error(log_lags('', 20))

expect_error(log_lags(20, ''))



expect_equal(log_lags(1.2, 4.4), log_lags(1, 4))

expect_equal(log_lags(10, 5), 0L:5L)

expect_equal(length(log_lags(15, 20)), 15L)

expect_equal(log_lags(15, 86400), sort(log_lags(15, 86400)))

expect_equal(log_lags(1, 86400), 0L)



expect_true(all(log_lags(2, 86400) <= 86400))

expect_true(all(log_lags(1, 86400) <= 86400))

expect_true(all(log_lags(86400, 86400) <= 86400))
