formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = rnorm(rows))


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# frecipes version
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepScale$new(x))$
  plate("tbl")
# recipes version
rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_scale(x) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

tinytest::expect_equivalent(frec, rec, info = "StepScale with R6 api")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# R6 version
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepScale$new(x))$
  plate("tbl")
# standard version
rec  = recipe(formula = formula, data = dat) |>
  step_scale(x) |>
  plate()

tinytest::expect_equivalent(frec, rec, info = "StepScale with recipes api")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

formula <- as.formula(y~x+z)
# frecipes version
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepScale$new(c(x,z)))$
  plate("tbl")
# recipes version
rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_scale(x, z) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

tinytest::expect_equivalent(frec, rec, info = "StepScale with multiple values")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
