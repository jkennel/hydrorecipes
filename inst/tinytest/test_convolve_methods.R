x <- rnorm(1000)
y <- rnorm(30)
y1 <- rev(y)

a <- frecipes:::convolve_filter(x, y, TRUE, TRUE)
b <- frecipes:::convolve_overlap_add(x, y1)
c <- frecipes:::convolve_overlap_save(x, y, 0)



tinytest::expect_equivalent(a, b)
tinytest::expect_equivalent(b, c)

x <- rnorm(1000 + 1)

a <- frecipes:::convolve_filter(x, y, TRUE, TRUE)
b <- frecipes:::convolve_overlap_add(x, y1)
c <- frecipes:::convolve_overlap_save(x, y, 0)


tinytest::expect_equivalent(a, b)
tinytest::expect_equivalent(b, c)

x <- rnorm(1000 + 1)
y <- rnorm(30 + 1)
y1 <- rev(y)

a <- frecipes:::convolve_filter(x, y, TRUE, TRUE)
b <- frecipes:::convolve_overlap_add(x, y1)
c <- frecipes:::convolve_overlap_save(x, y, 0)


tinytest::expect_equivalent(a, b)
tinytest::expect_equivalent(b, c)

x <- rnorm(1000)
y <- rnorm(30 + 1)
y1 <- rev(y)

a <- frecipes:::convolve_filter(x, y, TRUE, TRUE)
b <- frecipes:::convolve_overlap_add(x, y1)
c <- frecipes:::convolve_overlap_save(x, y, 0)


tinytest::expect_equivalent(a, b)
tinytest::expect_equivalent(b, c)

x <- rnorm(1000)
y <- rnorm(1001)
y1 <- rev(y)

tinytest::expect_error(frecipes:::convolve_filter(x, y, TRUE, TRUE))
tinytest::expect_error(frecipes:::convolve_overlap_add(x, y1))
tinytest::expect_error(frecipes:::convolve_overlap_save(x, y, 0))



# list versions -----------------------------------------------------------



n <- 864
x <- rnorm(n*100)
n_knots <- 9
max_lag <- 1 + n
knots <- frecipes:::log_lags_arma(n_knots, max_lag)
one_n <- c(1, length(knots))

y1 <- frecipes:::n_spline_list(0:n, 0L, 3L, knots[-one_n],
              knots[one_n], TRUE, FALSE,
              0L, FALSE)

bench::mark(
a <- frecipes:::convolve_list(x, y1, TRUE, TRUE),
c <- frecipes:::convolve_overlap_save_list(x, y1, 0)
)


tinytest::expect_equivalent(a, c)

