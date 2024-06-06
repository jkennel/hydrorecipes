
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
a <- hydrorecipes:::theis_aniso_time(distance_x = 10,
                                distance_y = 10,
                                storativity = 1e-6,
                                transmissivity_x = 1e-4,
                                transmissivity_y = 1e-4,
                                thickness = 1.0,
                                time = 1:10,
                                flow_rate = rep(1, 10))
b <- hydrorecipes:::grf_time(radius = sqrt(200),
                        specific_storage = 1e-6,
                        hydraulic_conductivity = 1e-4,
                        thickness = 1,
                        time = 1:10,
                        flow_rate = rep(1, 10),
                        flow_dimension = 2)

expect_equivalent(a, b, info = "Isotropic works")


dat <- data.frame(x = as.numeric(1:200),
                  y = rep(0.01, 200))
formula <- as.formula(y~x)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec1 = recipe(formula = formula, data = dat) |>
  step_aquifer_theis_aniso(time = x,
                     flow_rate = y) |>
  plate("dt")

frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferTheisAniso$new(time = x,
                                               flow_rate = y))$
  plate("dt")


expect_equivalent(frec1, frec2,
                  info = "R6 and hydrorecipes api are equivalent")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


  # a1 <- hydrorecipes:::theis_aniso_time(distance_x = 0,
  #                                      distance_y = 10,
  #                                      storativity = 1e-6,
  #                                      transmissivity_x = 1e-5,
  #                                      transmissivity_y = 1e-5,
  #                                      thickness = 1.0,
  #                                      time = 1:100,
  #                                      flow_rate = rep(0.01, 100))
  # a2 <- hydrorecipes:::theis_aniso_time(distance_x = 10,
  #                                      distance_y = 10,
  #                                      storativity = 1e-6,
  #                                      transmissivity_x = 1e-5,
  #                                      transmissivity_y = 1e-5,
  #                                      thickness = 1.0,
  #                                      time = 1:100,
  #                                      flow_rate = rep(0.01, 100))
  #
  # b1 <- hydrorecipes:::theis_aniso_time(distance_x = 0,
  #                                      distance_y = 10,
  #                                      storativity = 1e-6,
  #                                      transmissivity_x = 1e-4,
  #                                      transmissivity_y = 1e-4,
  #                                      thickness = 1.0,
  #                                      time = 1:100,
  #                                      flow_rate = rep(0.01, 100))
  # b2 <- hydrorecipes:::theis_aniso_time(distance_x = 10,
  #                                       distance_y = 10,
  #                                       storativity = 1e-6,
  #                                       transmissivity_x = 1e-4,
  #                                       transmissivity_y = 1e-4,
  #                                       thickness = 1.0,
  #                                       time = 1:100,
  #                                       flow_rate = rep(0.01, 100))
  # plot(a1[[1]]/a2[[1]], type = 'l', log = 'y')
  # abline(h = 1.2)
  # plot(b1[[1]]/b2[[1]], type = 'l', log = 'y')

  #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
