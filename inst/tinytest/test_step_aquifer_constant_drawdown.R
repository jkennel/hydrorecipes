
m <- matrix(
  c(1.00000000E-05,    1.74747218E+02,
    2.00000000E-05,    1.59975693E+02,
    5.00000000E-05,    1.43771934E+02,
    7.50000000E-05,    1.37572053E+02,
    1.00000000E-04,    1.33478472E+02,
    2.00000000E-04,    1.24525184E+02,
    5.00000000E-04,    1.14347107E+02,
    7.50000000E-04,    1.10345900E+02,
    1.00000000E-03,    1.07669736E+02,
    2.00000000E-03,    1.01717515E+02,
    5.00000000E-03,    9.47776445E+01,
    7.50000000E-03,    9.19964273E+01,
    1.00000000E-02,    9.01189516E+01,
    2.00000000E-02,    8.58921461E+01,
    5.00000000E-02,    8.08722811E+01,
    7.50000000E-02,    7.88319189E+01,
    1.00000000E-01,    7.74450895E+01,
    2.00000000E-01,    7.42944242E+01,
    5.00000000E-01,    7.05003520E+01,
    7.50000000E-01,    6.89416206E+01,
    1.00000000E+00,    6.78765619E+01,
    2.00000000E+00,    6.54399938E+01,
    5.00000000E+00,    6.24740172E+01,
    7.50000000E+00,    6.12452613E+01,
    1.00000000E+01,    6.04022180E+01,
    2.00000000E+01,    5.84628328E+01,
    5.00000000E+01,    5.60817145E+01,
    7.50000000E+01,    5.50886260E+01,
    1.00000000E+02,    5.44050045E+01),
  ncol = 2, byrow = TRUE)

times <- m[, 1]
s <- 10
Tr <- 10 # transmissivity of aquifer, m^2/d
S <- 1e-5 # storage coefficient of aquifer, -
rw <- 0.15 # radius of well, m
prec <- 1e-8
jl <- frecipes:::jacob_lohman_laplace(times, s, rw, Tr, S, prec, 12L)

tinytest::expect_equivalent(jl, m[, 2], tolerance = 1e-4, "c++ version")




formula <- formula(times~.)
dat <- data.frame(times = m[,1])


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepAquiferConstantDrawdown$new(time = times,
                                           drawdown = 10,
                                           thickness = 10,
                                           radius_well = 0.15,
                                           specific_storage = 1e-6,
                                           hydraulic_conductivity = 1,
                                           n_terms = 12L))$
  plate("tbl")

tinytest::expect_equivalent(frec[[2]], m[, 2],
                            tolerance = 1e-4,
                            info = "R6 api")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec = recipe(formula = formula, data = dat) |>
  step_aquifer_constant_drawdown(time = times,
                                 drawdown = 10,
                                 thickness = 10,
                                 radius_well = 0.15,
                                 specific_storage = 1e-6,
                                 hydraulic_conductivity = 1,
                                 n_terms = 12L) |>
  plate("tbl")


tinytest::expect_equivalent(frec[[2]], m[, 2],
                            tolerance = 1e-4,
                            info = "frecipes api")

