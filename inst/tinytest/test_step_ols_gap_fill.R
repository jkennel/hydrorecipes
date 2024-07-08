
set.seed(123)
n <- 100000
frm <- formula(x ~ y + z)


x <- cumsum(rnorm(n))
dat <- data.table(x = x, y = x, z = as.numeric(1:n))
dat[, x := x + c(rep(20, n/2), rep(0, n / 2))]
dat[, x := x + 3.0 * sin(z * 1 / n)]
tmp <- copy(dat$x)

# Set value to NA.  These values will be estimated.
dat[60000:70000, x := NA_real_]

dat <- unclass(dat)


h = recipe(formula = frm, data = dat) |>
  step_find_interval(z, vec = c(0, n/2, n)) |>
  step_intercept() |>
  step_spline_b(z, df = 4) |>
  step_drop_columns(z)

hrec = recipe(formula = frm, data = dat) |>
  step_ols_gap_fill(c(x, y, z), recipe = h) |>
  plate()


expect_equivalent(hrec$ols_gap_fill_x[60000:70000], tmp[60000:70000], tolerance = 0.1)

# plot(hrec$ols_gap_fill_x[60000:70000], type = 'l')
# points(tmp[60000:70000], type = 'l', col = 'red')
