
# semi-confined -----------------------------------------------------------


check_rojstaczer_fig_3 <- function(Q_W_ratio, transmissivity) {
  data("rojstaczer_1988b_fig_3")

  thickness_confining  <- 10
  thickness_vadose     <- 1
  diffusivity_confining <- 0.05

  storage_aquifer      <- 1e-4
  storage_confining    <- 1e-4

  attenuation          <- 1.0
  loading_efficiency   <- 0.5

  radius_well          <- 0.10
  diffusivity_vadose   <- 0.1

  vals <- rojstaczer_1988b_fig_3[Q_div_W == Q_W_ratio]

  frequency <- vals$W * transmissivity / (2.0 * pi * radius_well^2)


  roj_1988 <- areal_rojstaczer_semiconfined(frequency,
                                            radius_well,
                                            transmissivity,
                                            storage_confining,
                                            storage_aquifer,
                                            diffusivity_confining,
                                            diffusivity_vadose,
                                            thickness_confining,
                                            thickness_vadose,
                                            loading_efficiency,
                                            attenuation)
  wh_gain <- which(vals$variable == "gain")
  wh_phase <- which(vals$variable == "phase")

  expect_equal(Mod(roj_1988[wh_gain]), vals[wh_gain]$response, tolerance = 0.005)
  # had to add a 360 degree shift
  expect_equal(unwrap(Arg(roj_1988[wh_phase])) * 180 / pi, vals[wh_phase]$response, tolerance = 0.003)

}

Q_W_ratio <- 10000
transmissivity <- 1e-1
check_rojstaczer_fig_3(Q_W_ratio, transmissivity)

Q_W_ratio <- 1000
transmissivity <- 1e-2
check_rojstaczer_fig_3(Q_W_ratio, transmissivity)

Q_W_ratio <- 100
transmissivity <- 1e-3
check_rojstaczer_fig_3(Q_W_ratio, transmissivity)

Q_W_ratio <- 10
transmissivity <- 1e-4
check_rojstaczer_fig_3(Q_W_ratio, transmissivity)

Q_W_ratio <- 1
transmissivity <- 1e-5
check_rojstaczer_fig_3(Q_W_ratio, transmissivity)



check_rojstaczer_fig_5 <- function(storage_aquifer, storage_confining) {
  data("rojstaczer_1988b_fig_5")

  thickness_confining  <- 100
  thickness_vadose     <- 3
  diffusivity_confining <- 0.005

  transmissivity       <- 1e-1

  attenuation          <- 1.0
  loading_efficiency   <- 0.5

  radius_well          <- 0.10
  diffusivity_vadose   <- 0.1



  vals <- rojstaczer_1988b_fig_5[S == storage_confining]

  frequency <- vals$W * transmissivity / (2.0 * pi * radius_well^2)


  roj_1988 <- areal_rojstaczer_semiconfined(frequency,
                                            radius_well,
                                            transmissivity,
                                            storage_confining,
                                            storage_aquifer,
                                            diffusivity_confining,
                                            diffusivity_vadose,
                                            thickness_confining,
                                            thickness_vadose,
                                            loading_efficiency,
                                            attenuation)
  wh_gain <- which(vals$variable == "gain")
  wh_phase <- which(vals$variable == "phase")

  expect_equal(Mod(roj_1988[wh_gain]), vals[wh_gain]$response, tolerance = 0.005)
  expect_equal(unwrap(Arg(roj_1988[wh_phase])) * 180 / pi - 360, vals[wh_phase]$response, tolerance = 0.2)

}

storage_confining <- 1e-7
storage_aquifer  <- 1e-7
check_rojstaczer_fig_5(storage_aquifer, storage_confining)


storage_confining    <- 1e-5
storage_aquifer    <- 1e-5
check_rojstaczer_fig_5(storage_aquifer, storage_confining)

storage_confining    <- 1e-3
storage_aquifer    <- 1e-3
check_rojstaczer_fig_5(storage_aquifer, storage_confining)




check_rojstaczer_fig_5 <- function(R_div_Q, diffusivity_vadose) {

  data("rojstaczer_1988b_fig_6")

  thickness_confining  <- 10
  thickness_vadose     <- 100
  diffusivity_confining <- 0.001

  transmissivity       <- 1e-1

  attenuation          <- 1.0
  loading_efficiency   <- 0.5

  radius_well          <- 0.10
  storage_aquifer      <- 1e-7
  storage_confining    <- 1e-7


  vals <- rojstaczer_1988b_fig_6[R_div_Q == R_div_Q]
  frequency <- vals$dimensionless_frequency * 2.0 * diffusivity_confining / (2.0 * pi * thickness_confining^2)

  roj_1988 <- areal_rojstaczer_semiconfined(frequency,
                                            radius_well,
                                            transmissivity,
                                            storage_confining,
                                            storage_aquifer,
                                            diffusivity_confining,
                                            diffusivity_vadose,
                                            thickness_confining,
                                            thickness_vadose,
                                            loading_efficiency,
                                            attenuation)

  wh_gain <- which(vals$variable == "gain")
  wh_phase <- which(vals$variable == "phase")

  expect_equal(Mod(roj_1988[wh_gain]), vals[wh_gain]$response, tolerance = 0.005)
  expect_equal(unwrap(Arg(roj_1988[wh_phase])) * 180 / pi, vals[wh_phase]$response, tolerance = 0.2)

}

diffusivity_vadose   <- 0.0001
R_div_Q <- 1000
check_rojstaczer_fig_5(R_div_Q, diffusivity_vadose)

diffusivity_vadose   <- 0.001
R_div_Q <- 100
check_rojstaczer_fig_5(R_div_Q, diffusivity_vadose)

diffusivity_vadose   <- 0.01
R_div_Q <- 10
check_rojstaczer_fig_5(R_div_Q, diffusivity_vadose)

diffusivity_vadose   <- 0.1
R_div_Q <- 1
check_rojstaczer_fig_5(R_div_Q, diffusivity_vadose)

diffusivity_vadose   <- 100
R_div_Q <- 0.00011
check_rojstaczer_fig_5(R_div_Q, diffusivity_vadose)



data("rojstaczer_1988b_fig_6")

thickness_confining  <- 10
thickness_vadose     <- 100
diffusivity_confining <- 0.001

transmissivity       <- 1e-1

attenuation          <- 1.0
loading_efficiency   <- 0.5

radius_well          <- 0.10
storage_aquifer      <- 1e-7
storage_confining    <- 1e-7

diffusivity_vadose   <- 0.0001
R_div_Q <- 1000

vals <- rojstaczer_1988b_fig_6[R_div_Q == R_div_Q]
frequency <- vals$dimensionless_frequency * 2.0 * diffusivity_confining / (2.0 * pi * thickness_confining^2)

roj_1988_0 <- areal_rojstaczer_semiconfined(frequency,
                                            radius_well,
                                            transmissivity,
                                            storage_confining,
                                            storage_aquifer,
                                            diffusivity_confining,
                                            diffusivity_vadose,
                                            thickness_confining,
                                            thickness_vadose,
                                            loading_efficiency,
                                            attenuation)
formula <- as.formula(frequency~.)

# test api
frec1 = hydrorecipes:::Recipe$new(formula = formula, data = list(frequency = frequency))$
  add_step(hydrorecipes:::StepBaroFrequencySemiConfined$new(frequency,
                                                            radius_well,
                                                            transmissivity,
                                                            storage_confining,
                                                            storage_aquifer,
                                                            diffusivity_confining,
                                                            diffusivity_vadose,
                                                            thickness_confining,
                                                            thickness_vadose,
                                                            loading_efficiency,
                                                            attenuation))$
  plate("dt")

frec2 = recipe(formula = formula, data = list(frequency = frequency)) |>
  step_baro_frequency_semi_confined(frequency,
                                    radius_well,
                                    transmissivity,
                                    storage_confining,
                                    storage_aquifer,
                                    diffusivity_confining,
                                    diffusivity_vadose,
                                    thickness_confining,
                                    thickness_vadose,
                                    loading_efficiency,
                                    attenuation) |>
  plate("dt")
expect_equivalent(frec1[[2]], frec2[, 2],
                  info = "R6 and hydrorecipes api are equivalent")

expect_equivalent(frec1[[2]], roj_1988_0)





