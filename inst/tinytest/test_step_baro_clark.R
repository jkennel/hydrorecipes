set.seed(123)
n <- 1000
baro <- sin(seq(0, 10 * pi, length.out = n))
wl <- -0.4 * baro + rnorm(n, sd = 0.01)
clark <- be_clark_cpp(wl, baro, lag_space = 1, inverse = TRUE)

tinytest::expect_equivalent(clark, 0.40,  tolerance = 1e-2,
                            info = "be_clark_cpp works")




