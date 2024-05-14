formula <- as.formula(y~x)
rows <- 1000

dat <- data.frame(x = rnorm(rows),
                  y = collapse::qF(sample(1:10, rows, replace = TRUE)))

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# hydrorecipes version
frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepDummy$new(y))$
  plate("tbl")
# recipes version
rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_dummy(y, keep_original_cols = TRUE, one_hot = FALSE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)
expect_equivalent(frec, rec)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# hydrorecipes version
frec = recipe(formula = formula, data = dat) |>
  step_dummy(y, one_hot = FALSE) |>
  plate("tbl")
# recipes version
rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_dummy(y, keep_original_cols = TRUE, one_hot = FALSE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)
expect_equivalent(frec, rec)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# one hot
frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepDummy$new(y, one_hot = TRUE))$
  plate("tbl")
rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_dummy(y, keep_original_cols = TRUE, one_hot = TRUE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)
expect_equivalent(frec, rec)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
