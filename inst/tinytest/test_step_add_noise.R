formula <- as.formula(y~z)
rows <- 100000
set.seed(123)

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = rnorm(rows),
                  w = rnorm(rows))


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# hydrorecipes version
frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAddNoise$new(z))$
  plate("tbl")

expect_equivalent(sd(frec$z - dat$z), 1.0, tolerance = 0.01)
expect_equivalent(mean(frec$z - dat$z), 0.0, tolerance = 0.001)

frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAddNoise$new(c(y, z)))$
  plate("tbl")

expect_equivalent(sd(frec$y - dat$y), 1.0, tolerance = 0.01)
expect_equivalent(sd(frec$z - dat$z), 1.0, tolerance = 0.01)
expect_equivalent(mean(frec$y - dat$y), 0.0, tolerance = 0.01)
expect_equivalent(mean(frec$z - dat$z), 0.0, tolerance = 0.01)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
set.seed(123)
frec_1 = recipe(formula = formula, data = dat) |>
  step_add_noise(z) |>
  plate("tbl")
set.seed(123)
frec_2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAddNoise$new(z))$
  plate("tbl")

expect_equivalent(frec_1, frec_2)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

