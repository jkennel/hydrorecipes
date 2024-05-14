formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = as.numeric(1:rows),
                  z = rnorm(rows))

frec1 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepDropColumns$new(x))$
  plate("df")

frec2 = recipe(formula = formula, data = dat) |>
  step_drop_columns(x) |>
  plate("df")

expect_equivalent(frec1, frec2,
                            info = "R6 and hydrorecipes api are equivalent")

expect_equivalent(ncol(frec1), 1,
                            info = "R6 and hydrorecipes api are equivalent")
