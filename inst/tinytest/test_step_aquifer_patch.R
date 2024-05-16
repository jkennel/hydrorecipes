dat <- data.frame(x = as.numeric(1:200),
                  y = rep(0.01, 200))
formula <- as.formula(y~x)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec1 = recipe(formula = formula, data = dat) |>
  step_aquifer_patch(time = x,
                     flow_rate = 0.01,
                     thickness = 1.0,
                     radius = 100.0,
                     radius_patch = 200.0,
                     specific_storage_inner = 1e-6,
                     specific_storage_outer = 1e-6,
                     hydraulic_conductivity_inner = 1e-4,
                     hydraulic_conductivity_outer = 1e-4,
                     n_stehfest = 8L) |>
  plate("dt")

frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferPatch$new(time = x,
                                           flow_rate = 0.01,
                                           thickness = 1.0,
                                           radius = 100.0,
                                           radius_patch = 200.0,
                                           specific_storage_inner = 1e-6,
                                           specific_storage_outer = 1e-6,
                                           hydraulic_conductivity_inner = 1e-4,
                                           hydraulic_conductivity_outer = 1e-4,
                                           n_stehfest = 8L))$
  plate("dt")


expect_equivalent(frec1, frec2,
                            info = "R6 and hydrorecipes api are equivalent")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec1 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferPatch$new(time = x,
                                           flow_rate = 0.01,
                                           thickness = 1.0,
                                           radius = 100.0,
                                           radius_patch = 200.0,
                                           specific_storage_inner = 1e-6,
                                           specific_storage_outer = 1e-6,
                                           hydraulic_conductivity_inner = 1e-4,
                                           hydraulic_conductivity_outer = 1e-4,
                                           n_stehfest = 12L))$
  plate("dt")

frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferTheis$new(time = x,
                                           flow_rate = y,
                                           thickness = 1.0,
                                           radius = 100.0,
                                           specific_storage = 1.0e-6,
                                           hydraulic_conductivity = 1.0e-4))$
  plate("dt")

expect_equivalent(frec1, frec2,
                  info = "StepAquiferPatch is equivalent to StepAquiferTheis with inner and outer zones having the same properties", tolerance = 1e-4)


frec1 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferPatch$new(time = x,
                                               flow_rate = 0.01,
                                               thickness = 1.0,
                                               radius = 300.0,
                                               radius_patch = 200.0,
                                               specific_storage_inner = 1e-6,
                                               specific_storage_outer = 1e-6,
                                               hydraulic_conductivity_inner = 1e-4,
                                               hydraulic_conductivity_outer = 1e-4,
                                               n_stehfest = 12L))$
  plate("dt")

frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferTheis$new(time = x,
                                               flow_rate = y,
                                               thickness = 1.0,
                                               radius = 300.0,
                                               specific_storage = 1.0e-6,
                                               hydraulic_conductivity = 1.0e-4))$
  plate("dt")

# plot(x = frec1[[1]], y = frec1[[3]])
# points(x = frec2[[1]], y = frec2[[3]], type = 'l', col = 'red')

expect_equivalent(frec1, frec2,
                            info = "StepAquiferPatch is equivalent to StepAquiferTheis with inner and outer zones having the same properties", tolerance = 1e-4)
