# Inputs based on Chris Neville's work

m <- matrix(
  c(1.50000000E-01,1.00000000E-05,1.41143286E-03,
    1.50000000E-01,2.00000000E-05,2.81775134E-03,
    1.50000000E-01,5.00000000E-05,7.01133547E-03,
    1.50000000E-01,7.50000000E-05,1.04797309E-02,
    1.50000000E-01,1.00000000E-04,1.39259590E-02,
    1.50000000E-01,2.00000000E-04,2.75061962E-02,
    1.50000000E-01,5.00000000E-04,6.65127360E-02,
    1.50000000E-01,7.50000000E-04,9.72603994E-02,
    1.50000000E-01,1.00000000E-03,1.26576161E-01,
    1.50000000E-01,2.00000000E-03,2.31538441E-01,
    1.50000000E-01,5.00000000E-03,4.60951113E-01,
    1.50000000E-01,7.50000000E-03,5.88225503E-01,
    1.50000000E-01,1.00000000E-02,6.79495765E-01,
    1.50000000E-01,2.00000000E-02,8.67614707E-01,
    1.50000000E-01,5.00000000E-02,1.01411373E+00,
    1.50000000E-01,7.50000000E-02,1.05763794E+00,
    1.50000000E-01,1.00000000E-01,1.08552946E+00,
    1.50000000E-01,2.00000000E-01,1.14767186E+00,
    1.50000000E-01,5.00000000E-01,1.22463559E+00,
    1.50000000E-01,7.50000000E-01,1.25781130E+00,
    1.50000000E-01,1.00000000E+00,1.28116571E+00,
    1.50000000E-01,2.00000000E+00,1.33703352E+00,
    1.50000000E-01,5.00000000E+00,1.41039527E+00,
    1.50000000E-01,7.50000000E+00,1.44276568E+00,
    1.50000000E-01,1.00000000E+01,1.46571294E+00,
    1.50000000E-01,2.00000000E+01,1.52095686E+00,
    1.50000000E-01,5.00000000E+01,1.59392999E+00,
    1.50000000E-01,7.50000000E+01,1.62621058E+00,
    1.50000000E-01,1.00000000E+02,1.64911172E+00
  ),
  ncol = 3, byrow = TRUE)


times <- m[, 2]
Tr = 10 # transmissivity of aquifer, m^2/d
S = 1e-4 # storage coefficient of aquifer, -
rw = 0.15 # radius of well, m
rc = 0.15
r = 0.15
Q = 10
prec <- 1e-8
pc <- hydrorecipes:::papadopulos_cooper_laplace(times,
                                           Q,
                                           r,
                                           rc,
                                           rw,
                                           Tr,
                                           S,
                                           prec,
                                           12L)

expect_equivalent(pc, m[, 3], tolerance = 5e-4)


dat <- data.frame(x = m[, 2])
formula <- as.formula(.~x)

frec = recipe(formula = formula, data = dat) |>
  step_aquifer_wellbore_storage(time = x,
                                flow_rate = 10.0,
                                hydraulic_conductivity = 10.0,
                                specific_storage = 1e-4) |>
  prep() |>
  bake()

expect_equivalent(frec$result$aquifer_wellbore_storage, m[, 3],
                  tolerance = 5e-4)


formula <- as.formula(flow_rate~time)
q <- 0.05
t_1 <- 10000.0
t_2 <- 1.0
s_1 <- 1
s_2 <- 1e-5
r <- 10
R <- 0.15
time <- 10^seq(-4, 4, 0.1)
dat <- data.frame(time = time)

bh = recipe(formula = formula, data = dat) |>
  step_aquifer_patch(time = time,
                     flow_rate = q,
                     thickness = 1.0,
                     radius = r,
                     radius_patch = R,
                     specific_storage_inner = s_1,
                     specific_storage_outer = s_2,
                     hydraulic_conductivity_inner = t_1,
                     hydraulic_conductivity_outer = t_2,
                     n_stehfest = 12L) |>
  plate()


pc = recipe(formula = formula, data = dat) |>
  step_aquifer_wellbore_storage(time = time,
                                radius_casing = R,
                                radius_well = R,
                                radius = r,
                                flow_rate = q,
                                hydraulic_conductivity = t_2,
                                specific_storage = s_2) |>
  plate()

expect_equivalent(pc$aquifer_wellbore_storage, bh$aquifer_patch,
                  tolerance = 1e-6,
                  info = "StepAquiferWellboreStorage is equivalent to StepAquiferPatch with inner and outer zones having the same properties",
)




