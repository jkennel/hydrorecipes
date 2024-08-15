# unconfined --------------------------------------------------------------


check_rojstaczer_1990_fig_4 <- function(R_Qu_rat, diffusivity_vadose) {

  data('rojstaczer_1990_fig_4')

  thickness_aquifer    <- 100
  thickness_vadose     <- thickness_aquifer
  thickness_saturated_well <- thickness_aquifer
  diffusivity_vertical <- 0.01

  specific_yield       <- 0.02
  storage_aquifer      <- 0.001
  k_vertical           <- 0.0     # set ohm = 0

  transmissivity       <- 1e-4

  attenuation          <- 1.0
  loading_efficiency   <- 0.5


  vals <- rojstaczer_1990_fig_4[R_div_Qu == R_Qu_rat]

  frequency <- vals$Qu * 2.0 * diffusivity_vertical / (2.0 * pi * thickness_saturated_well^2)
  roj_1990 <- areal_rojstaczer_unconfined(frequency,
                                          radius_well,
                                          storage_aquifer,
                                          specific_yield,
                                          k_vertical,
                                          diffusivity_vertical,
                                          diffusivity_vadose,
                                          thickness_saturated_well,
                                          thickness_vadose,
                                          thickness_aquifer,
                                          loading_efficiency,
                                          attenuation)

  wh_gain <- which(vals$variable == "gain")
  wh_phase <- which(vals$variable == "phase")

  expect_equal(Mod(roj_1990[wh_gain]), vals[wh_gain]$response, tolerance = 0.02)
  expect_equal(unwrap(Arg(roj_1990[wh_phase])) * 180 / pi, vals[wh_phase]$response, tolerance = 0.2)

  return(NULL)

}


R_Qu_ratio <- 1000
diffusivity_vadose   <- 0.00001
check_rojstaczer_1990_fig_4(R_Qu_ratio, diffusivity_vadose)

R_Qu_ratio <- 100
diffusivity_vadose   <- 0.0001
check_rojstaczer_1990_fig_4(R_Qu_ratio, diffusivity_vadose)

R_Qu_ratio <- 10
diffusivity_vadose   <- 0.001
check_rojstaczer_1990_fig_4(R_Qu_ratio, diffusivity_vadose)

R_Qu_ratio <- 1
diffusivity_vadose   <- 0.01
check_rojstaczer_1990_fig_4(R_Qu_ratio, diffusivity_vadose)

R_Qu_ratio <- 1e-4
diffusivity_vadose   <- 10
check_rojstaczer_1990_fig_4(R_Qu_ratio, diffusivity_vadose)





data('rojstaczer_1990_fig_4')

thickness_aquifer    <- 100
thickness_vadose     <- thickness_aquifer
thickness_saturated_well <- thickness_aquifer
diffusivity_vertical <- 0.01

specific_yield       <- 0.02
storage_aquifer      <- 0.001
k_vertical           <- 0.0     # set ohm = 0

transmissivity       <- 1e-4
radius_well          <- 0.10

attenuation          <- 1.0
loading_efficiency   <- 0.5
R_Qu_ratio <- 1000
diffusivity_vadose   <- 0.00001

vals <- rojstaczer_1990_fig_4[R_div_Qu == R_Qu_ratio]
frequency <- vals$Qu * 2.0 * diffusivity_vertical / (2.0 * pi * thickness_saturated_well^2)

roj_1990_0 <- areal_rojstaczer_unconfined(frequency,
                                          radius_well,
                                          storage_aquifer,
                                          specific_yield,
                                          k_vertical,
                                          diffusivity_vertical,
                                          diffusivity_vadose,
                                          thickness_saturated_well,
                                          thickness_vadose,
                                          thickness_aquifer,
                                          loading_efficiency,
                                          attenuation)

formula <- as.formula(frequency~.)

# test api
frec1 = hydrorecipes:::Recipe$new(formula = formula, data = list(frequency = frequency))$
  add_step(hydrorecipes:::StepBaroFrequencyUnconfined$new(frequency,
                                                          radius_well,
                                                          storage_aquifer,
                                                          specific_yield,
                                                          k_vertical,
                                                          diffusivity_vertical,
                                                          diffusivity_vadose,
                                                          thickness_saturated_well,
                                                          thickness_vadose,
                                                          thickness_aquifer,
                                                          loading_efficiency,
                                                          attenuation))$
  plate("dt")

frec2 = recipe(formula = formula, data = list(frequency = frequency)) |>
  step_baro_frequency_unconfined(frequency,
                                 radius_well,
                                 storage_aquifer,
                                 specific_yield,
                                 k_vertical,
                                 diffusivity_vertical,
                                 diffusivity_vadose,
                                 thickness_saturated_well,
                                 thickness_vadose,
                                 thickness_aquifer,
                                 loading_efficiency,
                                 attenuation) |>
  plate("dt")
expect_equivalent(frec1[[2]], frec2[, 2],
                  info = "R6 and hydrorecipes api are equivalent")


expect_equivalent(frec1[[2]], roj_1990_0)

