
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




#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
dat <- data.frame(x = as.numeric(1:200),
                  y = rep(0.01, 200))
formula <- as.formula(y~x)
frec1 = recipe(formula = formula, data = dat) |>
  step_aquifer_theis_aniso(time = x,
                           flow_rate = y,
                           distance_x = 100,
                           distance_y = -100,
                           hydraulic_conductivity_major = 1e-5,
                           hydraulic_conductivity_minor = 1e-5,
                           major_axis_angle = 0) |>
  plate("dt")


frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepAquiferTheis$new(time = x,
                                               flow_rate = y,
                                               radius = sqrt(20000),
                                               hydraulic_conductivity = 1e-5))$
  plate("dt")


expect_equivalent(frec1, frec2,
                  info = "theis_aniso and theis give same results when anisotropy is 1")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


a <- theis_aniso_time(
  distance_x = 100,
  distance_y = 0.0,
  storativity = 1e-6,
  transmissivity_x = 1e-4,
  transmissivity_y = 1e-5,
  thickness = 1,
  time = dat$x,
  flow_rate = dat$y
)
b <- theis_aniso_time(
  distance_x = sqrt(5000),
  distance_y = sqrt(5000),
  storativity = 1e-6,
  transmissivity_x = 1e-4,
  transmissivity_y = 1e-5,
  thickness = 1,
  time = dat$x,
  flow_rate = dat$y
)

d <- theis_aniso_time(
  distance_x = 100,
  distance_y = 0.0,
  storativity = 1e-6,
  transmissivity_x = 1e-5,
  transmissivity_y = 1e-4,
  thickness = 1,
  time = dat$x,
  flow_rate = dat$y
)


freca = recipe(formula = formula, data = dat) |>
  step_aquifer_theis_aniso(time = x,
                           flow_rate = y,
                           distance_x = 100,
                           distance_y = 0,
                           hydraulic_conductivity_major = 1e-4,
                           hydraulic_conductivity_minor = 1e-5,
                           major_axis_angle = 0) |>
  plate("dt")

frecb = recipe(formula = formula, data = dat) |>
  step_aquifer_theis_aniso(time = x,
                           flow_rate = y,
                           distance_x = 100,
                           distance_y = 0,
                           hydraulic_conductivity_major = 1e-4,
                           hydraulic_conductivity_minor = 1e-5,
                           major_axis_angle = pi/4) |>
  plate("dt")
frecd = recipe(formula = formula, data = dat) |>
  step_aquifer_theis_aniso(time = x,
                           flow_rate = y,
                           distance_x = 100,
                           distance_y = 0,
                           hydraulic_conductivity_major = 1e-4,
                           hydraulic_conductivity_minor = 1e-5,
                           major_axis_angle = pi/2) |>
  plate("dt")






frecb1 = recipe(formula = formula, data = dat) |>
  step_aquifer_theis_aniso(time = x,
                           flow_rate = y,
                           distance_x = sqrt(5000),
                           distance_y = sqrt(5000),
                           hydraulic_conductivity_major = 1e-4,
                           hydraulic_conductivity_minor = 1e-5,
                           major_axis_angle = -pi/4) |>
  plate("dt")
frecb2 = recipe(formula = formula, data = dat) |>
  step_aquifer_theis_aniso(time = x,
                           flow_rate = y,
                           distance_x = sqrt(5000),
                           distance_y = sqrt(5000),
                           hydraulic_conductivity_major = 1e-4,
                           hydraulic_conductivity_minor = 1e-5,
                           major_axis_angle = pi/4) |>
  plate("dt")

b1 <- theis_aniso_time(
  distance_x = 100,
  distance_y = 0,
  storativity = 1e-6,
  transmissivity_x = 1e-4,
  transmissivity_y = 1e-5,
  thickness = 1,
  time = dat$x,
  flow_rate = dat$y
)
b2 <- theis_aniso_time(
  distance_x = 0,
  distance_y = 100,
  storativity = 1e-6,
  transmissivity_x = 1e-4,
  transmissivity_y = 1e-5,
  thickness = 1,
  time = dat$x,
  flow_rate = dat$y
)

expect_equivalent(a$generalized_radial, freca$aquifer_theis_aniso,
                  info = "test rotated versions")
expect_equivalent(b$generalized_radial, frecb$aquifer_theis_aniso,
                  info =  "test rotated versions")
expect_equivalent(d$generalized_radial, frecd$aquifer_theis_aniso,
                  tolerance = 1e-6,
                  info =  "test rotated versions")

expect_equivalent(b1$generalized_radial, frecb1$aquifer_theis_aniso,
                  tolerance = 1e-6,
                  info =  "test rotated versions")
expect_equivalent(b2$generalized_radial, frecb2$aquifer_theis_aniso,
                  tolerance = 1e-6,
                  info =  "test rotated versions")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
dat <- expand.grid(seq(0.1, 10, 0.01), seq(0.1, 10, 0.01))

# bench::mark(
#
# b1 <- hydrorecipes:::theis_aniso_grid(
#   distance_x = dat[,1],
#   distance_y = dat[,2],
#   storativity = 1e-6,
#   transmissivity_x = 1e-4,
#   transmissivity_y = 2e-5,
#   thickness = 1.0,
#   time = 100.0,
#   flow_rate = 0.01
# )
# )
#
# b1 <- matrix(b1$generalized_radial, ncol = 991)
# image(b1)
#
# b2 <- theis_aniso_time(
#   distance_x = 3,
#   distance_y = 3,
#   storativity = 1e-6,
#   transmissivity_x = 1e-4,
#   transmissivity_y = 1e-5,
#   thickness = 1,
#   time = 10,
#   flow_rate = 0.01
# )
