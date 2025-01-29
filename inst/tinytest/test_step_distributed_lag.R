formula <- as.formula(y~x)
rows <- 2e5

dat <- data.frame(x = rnorm(rows),
                  y = as.numeric(1:rows),
                  z = rnorm(rows))

frec1 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepDistributedLag$new(x,
                                             knots = hydrorecipes:::log_lags(6, 86401)))$
  prep()$
  bake()

frec2 = recipe(formula = formula, data = dat) |>
  step_distributed_lag(x, knots = hydrorecipes:::log_lags(6, 86401)) |>
  prep() |>
  bake()

expect_equivalent(frec1$result, frec2$result,
                  info = "R6 and hydrorecipes api are equivalent")


frec1$get_response_data(type = "df")


# spline checks
n <- 2e4
m <- sort(rnorm(n))
bk <- range(m)
knots <- collapse::fquantile(bk, probs = seq(0.05, 0.95, 0.3))

fr <- (collapse::qM(hydrorecipes:::b_spline_list(x = m,
                                             df = 0L,
                                             degree = 3L,
                                             internal_knots = knots,
                                             boundary_knots = bk,
                                             complete_basis = TRUE,
                                             periodic = FALSE,
                                             derivs = 0,
                                             integral = FALSE)))
bs <- unclass(splines2::bSpline(m,
                                knots = knots,
                                Boundary.knots = bk,
                                intercept = TRUE))

expect_equivalent(fr, bs,
                  info = "splines2 and hydrorecipes are equivalent (intercept)")

fr <- (collapse::qM(hydrorecipes:::b_spline_list(x = m,
                                   df = 0L,
                                   degree = 3L,
                                   internal_knots = knots,
                                   boundary_knots = bk,
                                   complete_basis = FALSE,
                                   periodic = FALSE,
                                   derivs = 0,
                                   integral = FALSE)))
bs <- unclass(splines2::bSpline(m,
                                knots = knots,
                                Boundary.knots = bk,
                                intercept = FALSE))

expect_equivalent(fr, bs,
                  info = "splines2 and hydrorecipes are equivalent (no intercept)")




fr <- (collapse::qM(hydrorecipes:::n_spline_list(x = m,
                                   df = 0L,
                                   degree = 3L,
                                   internal_knots = knots,
                                   boundary_knots = bk,
                                   complete_basis = FALSE,
                                   periodic = FALSE,
                                   derivs = 0,
                                   integral = FALSE)))
ns <- unclass(splines2::naturalSpline(m,
                                      knots = knots,
                                      Boundary.knots = bk,
                                      intercept = FALSE))
sns <- splines::ns(m,knots = knots, Boundary.knots = bk, intercept = FALSE)

expect_equivalent(fr, ns,
                  info = "splines2 and hydrorecipes are equivalent (no intercept)")


fr <- (collapse::qM(hydrorecipes:::n_spline_list(x = m,
                                   df = 0L,
                                   degree = 3L,
                                   internal_knots = knots,
                                   boundary_knots = bk,
                                   complete_basis = TRUE,
                                   periodic = FALSE,
                                   derivs = 0,
                                   integral = FALSE)))
ns <- unclass(splines2::naturalSpline(m,
                                      knots = knots,
                                      Boundary.knots = bk,
                                      intercept = TRUE))
sns <- splines::ns(m,knots = knots, Boundary.knots = bk, intercept = TRUE)

expect_equivalent(fr, ns,
                  info = "splines2 and hydrorecipes are equivalent (no intercept)")







# n <- 5e6
# m <- sort(rnorm(n))
# bk <- range(m)
# knots <- quantile(bk, probs = seq(0.05, 0.95, 0.1))
#
#
# bench::mark(
#
#   ns <- unclass(naturalSpline(m,
#                               knots = knots,
#                               Boundary.knots = bk,
#                               intercept = FALSE)),
#   bs <- unclass(bSpline(m,
#                         knots = knots,
#                         Boundary.knots = bk,
#                         intercept = FALSE)),
#   check = FALSE)
#
#
#
#
#
#
# set.seed(120)
# n <- 5
# m <- 0:6
# bk <- c(0,6)
# knots <- c(3,4,5)
#
# fr <- qM((hydrorecipes:::n_spline_list(x = m,
#                                  df = 0L,
#                                  degree = 3L,
#                                  internal_knots = knots,
#                                  boundary_knots = bk,
#                                  complete_basis = TRUE,
#                                  periodic = FALSE,
#                                  derivs = 0,
#                                  integral = FALSE)))
#
#
# n <- sort(rnorm(n * 5))
# mm <- qM(hydrorecipes::lag_list(n, 0:6, n_shift = 0, n_subset = 1))
#
# mm %*% (fr)
#
# formula = formula(n~.)
# dat <- list(n = n)
# frec2 = Recipe$new(formula = formula, data = dat)$
#   add_step(StepDistributedLag$new(n, knots = sort(c(bk, knots))))$
#   plate("df")
#
# frec1 = Recipe$new(formula = formula, data = dat)$
#   add_step(StepDistributedLag$new(n, basis_matrix = fr))$
#   plate("df")
#
# frec0 = Recipe$new(formula = formula, data = dat)$
#   add_step(StepDistributedLag$new(n, basis_matrix = qM(n_spline_list(0:6, 0L, 3L, knots,
#                                                                   bk, TRUE, FALSE,
#                                                                   0L, FALSE))))$
#   plate("df")
#
# rng = 0:self$lag_max
# one_n = c(1L, self$n_lag)
#
# n_spline_list(0:6, 0L, 3L, knots,
#               bk, TRUE, FALSE,
#               0L, FALSE)
#
# hydrorecipes:::convolve_list(n, fr, TRUE, TRUE)
#
