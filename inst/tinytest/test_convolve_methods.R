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


n <- 1e5
x <- rnorm(n)
y <- rnorm(n)
lag_max <- 1000

ccf_fft  <- na.omit(hydrorecipes:::convolve_correlation(as.numeric(x), as.numeric(y), lag_max))
ccf_base <- as.numeric(ccf(as.numeric(x), as.numeric(y), lag.max = lag_max, plot = FALSE)$acf)

expect_equivalent(ccf_fft,
                  ccf_base,
                  info = "stats::ccf and hydrorecipes:::convolve_correlation are the same (even, even)")


lag_max <- 1001

ccf_fft  <- na.omit(hydrorecipes:::convolve_correlation(as.numeric(x), as.numeric(y), lag_max))
ccf_base <- as.numeric(ccf(as.numeric(x), as.numeric(y), lag.max = lag_max, plot = FALSE)$acf)

expect_equivalent(ccf_fft, ccf_base,
                  info = "stats::ccf and hydrorecipes:::convolve_correlation are the same (even, odd)")


n <- 1e5 + 1
x <- rnorm(n)
y <- rnorm(n)
lag_max <- 1000

ccf_fft  <- na.omit(hydrorecipes:::convolve_correlation(as.numeric(x), as.numeric(y), lag_max))
ccf_base <- as.numeric(ccf(as.numeric(x), as.numeric(y), lag.max = lag_max, plot = FALSE)$acf)

expect_equivalent(ccf_fft, ccf_base,
                  info = "stats::ccf and hydrorecipes:::convolve_correlation are the same (odd, even)")


lag_max <- 1001

ccf_fft  <- na.omit(hydrorecipes:::convolve_correlation(as.numeric(x), as.numeric(y), lag_max))
ccf_base <- as.numeric(ccf(as.numeric(x), as.numeric(y), lag.max = lag_max, plot = FALSE)$acf)

expect_equivalent(ccf_fft, ccf_base,
                  info = "stats::ccf and hydrorecipes:::convolve_correlation are the same (odd, odd)")



# n <- 1e6 + 1
# x <- rnorm(n)
# y <- rnorm(n)
# lag_max <- 10000
#
# bench::mark(fft_ccf  <- hydrorecipes:::convolve_correlation(as.numeric(x), as.numeric(y), lag_max),
#             ccf_base <- as.numeric(ccf(as.numeric(x), as.numeric(y), lag.max = lag_max, plot = FALSE)$acf),
#             relative = FALSE
# )

# hydrorecipes version
dat <- data.frame(x, y)
form <- formula(x~y)
frec_1 = hydrorecipes:::Recipe$new(formula = form,
                                   data = dat)$
  add_step(hydrorecipes:::StepCrossCorrelation$new(c(x, y)))$
  prep()$
  bake()
frec_2 = recipe(formula = form, data = dat) |>
  step_cross_correlation(c(x, y)) |>
  prep() |>
  bake()

expect_equivalent(frec_1$steps[[2]]$acf,
                  frec_2$steps[[2]]$acf)

frec_1 = hydrorecipes:::Recipe$new(formula = form, data = dat)$
  add_step(hydrorecipes:::StepCrossCorrelation$new(x))$
  prep()$
  bake()
frec_2 = recipe(formula = form, data = dat) |>
  step_cross_correlation(x) |>
  prep() |>
  bake()

expect_equivalent(frec_1$steps[[2]]$acf,
                  frec_2$steps[[2]]$acf,
                  info = "R6 and standard method are the same")


