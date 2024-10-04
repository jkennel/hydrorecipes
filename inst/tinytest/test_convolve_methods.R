# even even
x <- rnorm(1000)
y <- rnorm(30)
y1 <- rev(y)

a <- hydrorecipes:::convolve_filter(x, y, TRUE, TRUE)
b <- hydrorecipes:::convolve_overlap_add(x, y1)
c <- hydrorecipes:::convolve_overlap_save(x, y, 0)

# even odd
x <- rnorm(1000)
y <- rnorm(30+1)
y1 <- rev(y)

a <- hydrorecipes:::convolve_filter(x, y, TRUE, TRUE)
b <- hydrorecipes:::convolve_overlap_add(x, y1)
c <- hydrorecipes:::convolve_overlap_save(x, y, 0)


expect_equivalent(a, b)
expect_equivalent(b, c)

# odd even
x <- rnorm(1000 + 1)
y <- rnorm(30)
y1 <- rev(y)

a <- hydrorecipes:::convolve_filter(x, y, TRUE, TRUE)
b <- hydrorecipes:::convolve_overlap_add(x, y1)
c <- hydrorecipes:::convolve_overlap_save(x, y, 0)


expect_equivalent(a, b)
expect_equivalent(b, c)

# odd odd
x <- rnorm(1000 + 1)
y <- rnorm(30 + 1)
y1 <- rev(y)

a <- hydrorecipes:::convolve_filter(x, y, TRUE, TRUE)
b <- hydrorecipes:::convolve_overlap_add(x, y1)
c <- hydrorecipes:::convolve_overlap_save(x, y, 0)


expect_equivalent(a, b)
expect_equivalent(b, c)


x <- rnorm(1000)
y <- rnorm(1001)
y1 <- rev(y)

expect_error(hydrorecipes:::convolve_filter(x, y, TRUE, TRUE))
expect_error(hydrorecipes:::convolve_overlap_add(x, y1))
expect_error(hydrorecipes:::convolve_overlap_save(x, y, 0))



# list versions -----------------------------------------------------------



n <- 864
x <- rnorm(n*100)
n_knots <- 9
max_lag <- 1 + n
knots <- hydrorecipes:::log_lags(n_knots, max_lag)
one_n <- c(1, length(knots))

y1 <- hydrorecipes:::n_spline_list(0:n, 0L, 3L, knots[-one_n],
              knots[one_n], TRUE, FALSE,
              0L, FALSE)

bench::mark(
a <- hydrorecipes:::convolve_list(x, y1, TRUE, TRUE),
c <- hydrorecipes:::convolve_overlap_save_list(x, y1, 0)
)


expect_equivalent(a, c)



# even even
x <- rnorm(10)
y <- rep(1/3, 3)
y1 <- rev(y)

c <- hydrorecipes:::convolve_overlap_save(x, y, 1)
d <- caTools::runmean(x, 3, align = "center", endrule = "NA")
expect_equivalent(c, d)



c <- hydrorecipes:::convolve_overlap_save(x, y, 0)
d <- caTools::runmean(x, 3, align = "right", endrule = "NA")
expect_equivalent(c, d)


c <- hydrorecipes:::convolve_overlap_save(x, y, 2)
d <- caTools::runmean(x, 3, align = "left", endrule = "NA")
expect_equivalent(c, d)

