set.seed(123)
n <- 2e6
baro <- sin(seq(0, 10 * pi, length.out = n))
wl <- -0.4 * baro + rnorm(n, sd = 0.01)
clark <- be_clark_cpp(wl, baro, lag_space = 1L, inverse = TRUE)

tinytest::expect_equivalent(clark, 0.40,  tolerance = 1e-2,
                            info = "be_clark_cpp works")




# frecipes version
dat <- data.frame(x = baro,
                  y = wl)
formula <- formula(y~x)
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepBaroClark$new(y,
                             x,
                             1L,
                             TRUE))$
  prep()$
  bake()


baro2 <- sin(seq(0, 10 * pi, length.out = n))
wl2 <- 0.4 * baro + rnorm(n, sd = 0.01)
dat2 <- data.frame(x = baro2,
                   y = wl2)
frec2 = Recipe$new(formula = formula, data = dat2)$
  add_step(StepBaroClark$new(y,
                             x,
                             1L,
                             FALSE))$
  prep()$
  bake()

tinytest::expect_equivalent(frec$steps[[1]]$barometric_efficiency,
                            frec2$steps[[1]]$barometric_efficiency,
                            tolerance = 1e-2,
                            info = "step_baro_clark inverse works")


frec2 = recipe(formula = formula, data = dat) |>
  step_baro_clark(y, x, 1L, TRUE) |>
  prep() |>
  bake()

tinytest::expect_equivalent(frec$steps[[1]]$barometric_efficiency,
                            frec2$steps[[1]]$barometric_efficiency,
                            tolerance = 1e-2,
                            info = "step_baro_clark R6 and frecipes api are equivalent")


frec2 = Recipe$new(formula = formula, data = dat2)$
  add_step(StepBaroClark$new(y,
                             x,
                             1L,
                             FALSE))$
  bake()


dat <- data.frame(x = baro,
                  y = wl)
formula <- formula(y~x)
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepBaroClark$new(y,
                             x,
                             1:5,
                             TRUE))$
  prep()$
  bake()

tinytest::expect_equivalent(frec$steps[[2]]$barometric_efficiency,
                            rep(0.4, 5),
                            tolerance = 1e-2,
                            info = "step_baro_clark multiple lag_space works")
