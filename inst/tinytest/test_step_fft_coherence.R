# formula <- as.formula(y~x)
#
# dat <- data.frame(x = rnorm(200),
#                   y = rnorm(200))
#
# formula <- as.formula(.~x+y)
# frec1 = Recipe$new(formula = formula, data = dat)$
#   add_step(StepPgram$new(c(x, y)))$
#   add_step(StepCoherence$new())$
#   plate("df")
#
# frec2 = recipe(formula = formula, data = dat) |>
#   step_fft_pgram(c(x, y)) |>
#   step_fft_coherence() |>
#   plate("df")
#
# expect_equivalent(frec1, frec2,
#                             info = "R6 and hydrorecipes api are equivalent")

