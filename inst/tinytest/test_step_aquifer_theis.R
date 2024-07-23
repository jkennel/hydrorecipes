dat <- data.frame(x = as.numeric(1:200),
                  y = rep(0.01, 200))
formula <- as.formula(y~x)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec1 = recipe(formula = formula, data = dat) |>
  step_aquifer_theis(time = x,
                     flow_rate = y) |>
  plate("dt")

frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferTheis$new(time = x,
                                flow_rate = y))$
  plate("dt")


expect_equivalent(frec1, frec2,
                            info = "R6 and hydrorecipes api are equivalent")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

dat <- data.frame(x = as.numeric(1:200),
                  y = rep(0.01, 200))
formula <- as.formula(y~x)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec1 = recipe(formula = formula, data = dat) |>
  step_aquifer_theis(time = x,
                     flow_rate = y,
                     hydraulic_conductivity = 1e-4,
                     specific_storage = 1e-6) |>
  plate("dt")
frec2 = recipe(formula = formula, data = dat) |>
  step_aquifer_grf(time = x,
                   flow_rate = y,
                   hydraulic_conductivity = 1e-4,
                   specific_storage = 1e-6) |>
  plate("dt")


expect_equivalent(frec1, frec2,
                  info = "GRF and Theis api are equivalent")


frec1 = recipe(formula = formula, data = dat) |>
  step_aquifer_grf(time = x,
                     flow_rate = y,
                     hydraulic_conductivity = 1e-4,
                     specific_storage = 1e-6,
                     prefix = "test_name_1") |>
  step_aquifer_grf(time = x,
                     flow_rate = y,
                     hydraulic_conductivity = 1e-4,
                     specific_storage = 1e-6,
                     prefix = "test_name_2") |>
  plate("dt")

frec1 <- names(frec1)[3:4]

expect_equivalent(frec1, c("test_name_1", "test_name_2"),
                  info = "Theis column naming works")

