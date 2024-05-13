formula <- as.formula(y~z)
rows <- 100000

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = rnorm(rows),
                  w = rnorm(rows))


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# frecipes version
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepAddNoise$new(z))$
  plate("tbl")

tinytest::expect_equivalent(sd(frec$z - dat$z), 1.0, tolerance = 0.01)
tinytest::expect_equivalent(mean(frec$z - dat$z), 0.0, tolerance = 0.001)

frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepAddNoise$new(c(y, z)))$
  plate("tbl")

tinytest::expect_equivalent(sd(frec$y - dat$y), 1.0, tolerance = 0.01)
tinytest::expect_equivalent(sd(frec$z - dat$z), 1.0, tolerance = 0.01)
tinytest::expect_equivalent(mean(frec$y - dat$y), 0.0, tolerance = 0.01)
tinytest::expect_equivalent(mean(frec$z - dat$z), 0.0, tolerance = 0.01)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
set.seed(123)
frec_1 = recipe(formula = formula, data = dat) |>
  step_add_noise(z) |>
  plate("tbl")
set.seed(123)
frec_2 = Recipe$new(formula = formula, data = dat)$
  add_step(StepAddNoise$new(z))$
  plate("tbl")

tinytest::expect_equivalent(frec_1, frec_2)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

