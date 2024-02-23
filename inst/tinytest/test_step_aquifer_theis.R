dat <- data.frame(x = as.numeric(1:200),
                  y = rep(0.01, 200))
formula <- as.formula(y~x)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec1 = recipe(formula = formula, data = dat) |>
  step_aquifer_theis(time = x,
                     flow_rate = y) |>
  plate("dt")

frec2 = Recipe$new(formula = formula, data = dat)$
  add_step(StepAquiferTheis$new(time = x,
                                flow_rate = y))$
  plate("dt")


tinytest::expect_equivalent(frec1, frec2,
                            info = "R6 and frecipes api are equivalent")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
