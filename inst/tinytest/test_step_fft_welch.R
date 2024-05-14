formula <- as.formula(y~.)

dat <- data.frame(x = rnorm(200),
                  y = rnorm(200))

formula <- as.formula(.~x+y)
frec1 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepWelch$new(c(x, y), length_subset = 10, window = hydrorecipes:::window_rectangle(10)))$
  plate("df")

frec2 = recipe(formula = formula, data = dat) |>
  step_fft_welch(c(x,y), length_subset = 10, window = hydrorecipes:::window_rectangle(10)) |>
  plate("df")

expect_equivalent(frec1, frec2,
                            info = "R6 and hydrorecipes api are equivalent")

