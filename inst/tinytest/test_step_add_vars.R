formula <- as.formula(y~z)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = rnorm(rows),
                  w = rnorm(rows))


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# hydrorecipes version
frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAddVars$new(c(x)))$
  plate("tbl")

expect_equivalent(ncol(frec), 3)

frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAddVars$new("x"))$
  plate("tbl")

expect_equivalent(ncol(frec), 3)


frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAddVars$new(c(x, w)))$
  plate("tbl")

expect_equivalent(ncol(frec), 4)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec = recipe(formula = formula, data = dat) |>
  step_add_vars(c(x, w)) |>
  plate("tbl")

expect_equivalent(ncol(frec), 4)

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

