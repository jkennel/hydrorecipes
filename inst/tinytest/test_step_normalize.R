formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = rnorm(rows))
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepNormalize$new(x))$
  prep()$
  bake()$
  data("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_normalize(x) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

tinytest::expect_equivalent(frec, rec)
