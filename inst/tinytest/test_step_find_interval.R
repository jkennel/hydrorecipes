formula <- as.formula(y~x)

rows = 200
dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = rnorm(rows))

frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepFindInterval$new(x, vec = c(-0.1, 0.0, 0.1)))$
  plate("tbl")

frec1 = recipe(formula = formula, data = dat) |>
  step_find_interval(x, vec = c(-0.1, 0.0, 0.1)) |>
  plate("tbl")

expect_equivalent(frec, frec1, info = "StepFindInterval with R6 api")
