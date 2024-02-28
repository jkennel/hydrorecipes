set.seed(123)
n <- 1000
baro <- sin(seq(0, 10 * pi, length.out = n))
wl <- -0.4 * baro + rnorm(n, sd = 0.01)

lsq_diff <- be_least_squares_diff_cpp(wl, baro, lag_space = 1, inverse = TRUE)
lsq <- be_least_squares_cpp(wl, baro, inverse = TRUE)

tinytest::expect_equivalent(lsq_diff, 0.40,  tolerance = 1e-3,
                            info = "be_least_squares_diff_cpp works")



tinytest::expect_equivalent(lsq, 0.40,  tolerance = 1e-3,
                            info = "be_least_squares_cpp works")


lsq_diff_1 <- be_least_squares_diff_cpp(wl, baro, lag_space = 1, inverse = TRUE)
lsq_diff_5 <- be_least_squares_diff_cpp(wl, baro, lag_space = 5, inverse = TRUE)

tinytest::expect_equivalent(lsq_diff_1, lsq_diff_5, tolerance = 1e-2,
                            info = "be_least_squares_diff_cpp works")

