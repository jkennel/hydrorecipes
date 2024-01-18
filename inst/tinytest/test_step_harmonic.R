formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = as.numeric(1:rows),
                  z = rnorm(rows))
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepHarmonic$new(y,
                            frequency = c(3),
                            cycle_size = 0.1,
                            starting_value = 0))$
  prep()$
  bake()$
  data("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_harmonic(y,
                         frequency = c(3),
                         cycle_size = 0.1,
                         starting_value = 0,
                         keep_original_cols = TRUE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

tinytest::expect_equivalent(frec, rec)
