dat <- data.frame(x = as.numeric(1:20),
                  y = rep(0.01, 20))
formula <- as.formula(y~x)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec1 = recipe(formula = formula, data = dat) |>
  step_aquifer_grf(time = x,
                   flow_rate = y) |>
  plate("dt")

frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferGRF$new(time = x,
                                             flow_rate = y))$
  plate("dt")


expect_equivalent(frec1, frec2,
                            info = "R6 and hydrorecipes api are equivalent")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec1 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferTheis$new(time = x,
                              flow_rate = y))$
  plate("dt")

frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferGRF$new(time = x,
                              flow_rate = y))$
  plate("dt")

expect_equivalent(frec1[[2]], frec2[, 2],
                            info = "Theis and GRF (radial) are equivalent")
