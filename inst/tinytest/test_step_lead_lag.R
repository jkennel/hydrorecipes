formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = as.numeric(1:rows),
                  z = rnorm(rows))
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepLeadLag$new(y, lag = 1))$
  plate("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_lag(y, lag = 1, keep_original_cols = TRUE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

tinytest::expect_equivalent(frec, rec)
