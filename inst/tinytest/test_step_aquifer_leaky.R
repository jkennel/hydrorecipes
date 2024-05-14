dat <- data.frame(x = as.numeric(1:200),
                  y = rep(0.01, 200))
formula <- as.formula(y~x)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec1 = recipe(formula = formula, data = dat) |>
  step_aquifer_leaky(time = x,
                     flow_rate = y) |>
  plate("dt")

frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferLeaky$new(time = x,
                                           flow_rate = y))$
  plate("dt")


expect_equivalent(frec1, frec2,
                            info = "R6 and hydrorecipes api are equivalent")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec1 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferLeaky$new(time = x,
                                           flow_rate = y,
                                           leakage = 10000000))$
  plate("dt")
frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferTheis$new(time = x,
                                           flow_rate = y))$
  plate("dt")

expect_equivalent(frec1, frec2,
                            info = "StepAquiferLeaky is equivalent to StepAquiferTheis when leakage is large")
