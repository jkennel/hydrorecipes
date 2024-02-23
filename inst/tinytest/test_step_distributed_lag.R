formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = as.numeric(1:rows),
                  z = rnorm(rows))

frec1 = Recipe$new(formula = formula, data = dat)$
  add_step(StepDistributedLag$new(x,
                                  knots = frecipes:::log_lags_arma(5, 86401)))$
  plate("df")

frec2 = recipe(formula = formula, data = dat) |>
  step_distributed_lag(x, knots = frecipes:::log_lags_arma(5, 86401)) |>
  plate("df")

tinytest::expect_equivalent(frec1, frec2,
                            info = "R6 and frecipes api are equivalent")
