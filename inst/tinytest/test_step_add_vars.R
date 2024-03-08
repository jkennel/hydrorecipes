formula <- as.formula(y~z)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = rnorm(rows),
                  w = rnorm(rows))


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# frecipes version
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepAddVars$new(x))$
  plate("tbl")

tinytest::expect_equivalent(ncol(frec),3)

frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepAddVars$new(c(x, w)))$
  plate("tbl")

tinytest::expect_equivalent(ncol(frec), 4)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec = recipe(formula = formula, data = dat) |>
  step_add_vars(x) |>
  plate("tbl")

tinytest::expect_equivalent(ncol(frec), 3)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

