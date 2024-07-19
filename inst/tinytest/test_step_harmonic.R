formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = as.numeric(1:rows),
                  z = rnorm(rows))
frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepHarmonic$new(y,
                                           frequency = c(3),
                                           cycle_size = 0.1,
                                           starting_value = 0))$
  plate("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_harmonic(y,
                         frequency = c(3),
                         cycle_size = 0.1,
                         starting_val = 0.0,
                         keep_original_cols = TRUE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

expect_equivalent(frec, rec[,c(1,2,4,3)])
