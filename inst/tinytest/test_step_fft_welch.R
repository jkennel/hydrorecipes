formula <- as.formula(y~.)

dat <- data.frame(x = rnorm(200),
                  y = rnorm(200))

formula <- as.formula(.~x+y)
frec1 = Recipe$new(formula = formula, data = dat)$
  add_step(StepWelch$new(c(x, y), length_subset = 10, window = window_rectangle(10)))$
  plate("df")

frec2 = recipe(formula = formula, data = dat) |>
  step_fft_welch(c(x,y), length_subset = 10, window = window_rectangle(10)) |>
  plate("df")

tinytest::expect_equivalent(frec1, frec2,
                            info = "R6 and frecipes api are equivalent")

