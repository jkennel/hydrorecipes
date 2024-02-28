# check vs. CBP 1967 table 1
Tt_r2    <- rep(c(1.0, 2.15, 4.64), 4) * 10^(rep(c(-3, -2, -1, 0), each = 3))
a_01     <- c(0.9771, 0.9658, 0.9490, 0.9238, 0.8860, 0.8283,
              0.7460, 0.6289, 0.4782, 0.3117, 0.1665, 0.07415)
a_0001   <- c(0.9969, 0.9949, 0.9914, 0.9853, 0.9744, 0.9545,
              0.9183, 0.8538, 0.7436, 0.5729, 0.3543, 0.1554)
a_000001 <- c(0.9992, 0.9985, 0.9970, 0.9942, 0.9888, 0.9781,
              0.9572, 0.9167, 0.8410, 0.7080, 0.5038, 0.2620)

dat <- list(x = Tt_r2)

formula = formula(x~.)
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

frec1 = recipe(formula = formula, data = dat) |>
  step_slug_cbp(
    times = x,
    radius = 1.0,
    radius_casing = 1.0,
    radius_well = 1.0,
    specific_storage = 1e-1,
    hydraulic_conductivity = 1.0,
    thickness = 1.0,
    head_0 = 1.0,
    n_terms = 12L
  ) |>
  plate("dt")




tinytest::expect_equivalent(frec1[[2]], a_01, tolerance = 1e-3,
                            info = "R6 and frecipes api are equivalent")

frec1 = recipe(formula = formula, data = dat) |>
  step_slug_cbp(
    times = x,
    radius = 1.0,
    radius_casing = 1.0,
    radius_well = 1.0,
    specific_storage = 1e-3,
    hydraulic_conductivity = 1.0,
    thickness = 1.0,
    head_0 = 1.0,
    n_terms = 12L
  ) |>
  plate("dt")

tinytest::expect_equivalent(frec1[[2]], a_0001, tolerance = 1e-3,
                            info = "R6 and frecipes api are equivalent")

frec1 = recipe(formula = formula, data = dat) |>
  step_slug_cbp(
    times = x,
    radius = 1.0,
    radius_casing = 1.0,
    radius_well = 1.0,
    specific_storage = 1e-5,
    hydraulic_conductivity = 1.0,
    thickness = 1.0,
    head_0 = 1.0,
    n_terms = 12L
  ) |>
  plate("dt")

tinytest::expect_equivalent(frec1[[2]], a_000001, tolerance = 1e-3,
                            info = "R6 and frecipes api are equivalent")

dat <- list(x = 0.0)
h_0 <- 20
frec1 = recipe(formula = formula, data = dat) |>
  step_slug_cbp(
    times = x,
    radius = 1.0,
    radius_casing = 1.0,
    radius_well = 1.0,
    specific_storage = 1e-5,
    hydraulic_conductivity = 1.0,
    thickness = 1.0,
    head_0 = h_0,
    n_terms = 12L
  ) |>
  plate("dt")

tinytest::expect_equivalent(frec1[[2]], h_0,
                            info = "R6 and frecipes api are equivalent")

