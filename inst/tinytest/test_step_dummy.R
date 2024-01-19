formula <- as.formula(y~x)
rows <- 1000

dat <- data.frame(x = rnorm(rows),
                  y = qF(sample(1:10, rows, replace = TRUE)))
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepDummy$new(y))$
  plate("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_dummy(y, keep_original_cols = TRUE, one_hot = FALSE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)
tinytest::expect_equivalent(frec, rec)


frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepDummy$new(y, one_hot = TRUE))$
  plate("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_dummy(y, keep_original_cols = TRUE, one_hot = TRUE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)
tinytest::expect_equivalent(frec, rec)
