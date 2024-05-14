# formula <- as.formula(y~z)
# rows <- 20
#
# dat <- data.frame(x = rnorm(rows),
#                   y = 1:rows,
#                   z = rnorm(rows),
#                   w = rnorm(rows))
#
#
# #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# # add noise ---------------------------------------------------------------
#
# frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
#   add_step(hydrorecipes:::StepAddNoise$new(x))$
#   prep()$
#   bake()
# frec$tidy()
# # add vars -----------------------------------------------------------------
#
# frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
#   add_step(hydrorecipes:::StepAddVars$new(x))$
#   prep()$
#   bake()
# frec$tidy()
#
# # constant drawdown-----------------------------------------------------------------
#
#
# m <- matrix(
#   c(1.00000000E-05,    1.74747218E+02,
#     2.00000000E-05,    1.59975693E+02,
#     5.00000000E-05,    1.43771934E+02,
#     7.50000000E-05,    1.37572053E+02,
#     1.00000000E-04,    1.33478472E+02,
#     2.00000000E-04,    1.24525184E+02,
#     5.00000000E-04,    1.14347107E+02,
#     7.50000000E-04,    1.10345900E+02,
#     1.00000000E-03,    1.07669736E+02,
#     2.00000000E-03,    1.01717515E+02,
#     5.00000000E-03,    9.47776445E+01,
#     7.50000000E-03,    9.19964273E+01,
#     1.00000000E-02,    9.01189516E+01,
#     2.00000000E-02,    8.58921461E+01,
#     5.00000000E-02,    8.08722811E+01,
#     7.50000000E-02,    7.88319189E+01,
#     1.00000000E-01,    7.74450895E+01,
#     2.00000000E-01,    7.42944242E+01,
#     5.00000000E-01,    7.05003520E+01,
#     7.50000000E-01,    6.89416206E+01,
#     1.00000000E+00,    6.78765619E+01,
#     2.00000000E+00,    6.54399938E+01,
#     5.00000000E+00,    6.24740172E+01,
#     7.50000000E+00,    6.12452613E+01,
#     1.00000000E+01,    6.04022180E+01,
#     2.00000000E+01,    5.84628328E+01,
#     5.00000000E+01,    5.60817145E+01,
#     7.50000000E+01,    5.50886260E+01,
#     1.00000000E+02,    5.44050045E+01),
#   ncol = 2, byrow = TRUE)
#
#
# formula <- formula(times~.)
# dat <- data.frame(times = m[,1])
#
# frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
#   add_step(hydrorecipes:::StepAquiferConstantDrawdown$new(time = times,
#                                            drawdown = 10,
#                                            thickness = 10,
#                                            radius_well = 0.15,
#                                            specific_storage = 1e-6,
#                                            hydraulic_conductivity = 1,
#                                            n_terms = 12L))$
#   prep()$
#   bake()
# frec$tidy()
#
#
# # grf ----------------------------------------------------------------
#
# dat <- data.frame(x = as.numeric(1:1e4),
#                   y = rep(0.01, 1e4))
#
# frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
#   add_step(hydrorecipes:::StepAquiferGRF$new(time = x,
#                               flow_rate = y))$
#   prep()$
#   bake()
# frec$tidy()
#
#
# # leaky --------------------------------------------------------------
# dat <- data.frame(x = as.numeric(1:1e4),
#                   y = rep(0.01, 1e4))
#
# frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
#   add_step(hydrorecipes:::StepAquiferLeaky$new(time = x,
#                                 flow_rate = y))$
#   prep()$
#   bake()
# frec$tidy()
#
# # aquifer patch--------------------------------------------------------------
#
# dat <- data.frame(x = as.numeric(1:200),
#                   y = rep(0.01, 200))
# formula <- as.formula(y~x)
#
#
# frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
#   add_step(hydrorecipes:::StepAquiferPatch$new(time = x,
#                                 flow_rate = 0.01,
#                                 thickness = 1.0,
#                                 radius = 100.0,
#                                 radius_patch = 200.0,
#                                 specific_storage_inner = 1e-6,
#                                 specific_storage_outer = 1e-6,
#                                 hydraulic_conductivity_inner = 1e-4,
#                                 hydraulic_conductivity_outer = 1e-4,
#                                 n_stehfest = 8L))$
#   prep()$
#   bake()
# frec$tidy()
#
# # theis ---------------------------------------------------------------
#
# dat <- data.frame(x = as.numeric(1:1e4),
#                   y = rep(0.01, 1e4))
#
# frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
#   add_step(hydrorecipes:::StepAquiferTheis$new(time = x,
#                               flow_rate = y))$
#   prep()$
#   bake()
# frec$tidy()
# # baro acworth --------------------------------------------------------------
# data("kennel_2020")
#
# # synthetic earthtide
# latitude     <- 34.23411                           # latitude
# longitude    <- -118.678                           # longitude
# elevation    <- 500                                # elevation
# cutoff       <- 1e-5                               # cutoff
# catalog      <- 'ksm04'                            # hartmann wenzel catalog
# astro_update <- 300                                # how often to update astro parameters
# method       <- 'volume_strain'                    # which potential to calculate
#
# wave_groups_dl <- as.data.table(earthtide::eterna_wavegroups)
# wave_groups_dl <- na.omit(wave_groups_dl[time == '1 month'])
# wave_groups_dl <- wave_groups_dl[wave_groups_dl$start > 0.5,]
# wave_groups_dl <- wave_groups_dl[, list(start, end)]
# ngr <- nrow(wave_groups_dl)
#
#
# frec = hydrorecipes:::Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
#   add_step(hydrorecipes:::StepEarthtide$new(datetime,
#                              wave_groups = wave_groups_dl,
#                              latitude = latitude,
#                              longitude = longitude,
#                              elevation = elevation,
#                              cutoff = cutoff,
#                              catalog = catalog))$
#   plate('dt')
#
# # frec[, wl := lm(wl~datetime)$residuals]
# # frec[, baro := lm(baro~datetime)$residuals]
# # frec[, earthtide := lm(earthtide~datetime)$residuals]
#
# # frec[, wl := wl - baro]
# # frec2 = hydrorecipes:::Recipe$new(formula = wl~baro+datetime+earthtide, data = frec)$
# #   add_step(hydrorecipes:::StepBaroAcworth$new(water_level = wl,
# #                                baro = baro,
# #                                earth_tide = earthtide,
# #                                spans = 3))$
# #   prep()$
# #   bake()
# #
# #
# # frec2 = hydrorecipes:::Recipe$new(formula = wl~baro+datetime+earthtide, data = frec)$
# #   add_step(hydrorecipes:::StepBaroRau$new(water_level = wl,
# #                                baro = baro,
# #                                earth_tide = earthtide,
# #                                time = datetime))$
# #   prep()$
# #   bake()
#
#
# # earthtide --------------------------------------------------------------
#
#
# frec = hydrorecipes:::Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
#   add_step(hydrorecipes:::StepEarthtide$new(datetime,
#                              wave_groups = wave_groups_dl,
#                              latitude = latitude,
#                              longitude = longitude,
#                              elevation = elevation,
#                              cutoff = cutoff,
#                              catalog = catalog))$
#   prep()$
#   bake()
# frec$tidy()
#
#
# frec = hydrorecipes:::Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
#   add_step(hydrorecipes:::StepEarthtide$new(datetime,
#                              do_predict = FALSE,
#                              wave_groups = wave_groups_dl,
#                              latitude = latitude,
#                              longitude = longitude,
#                              elevation = elevation,
#                              cutoff = cutoff,
#                              catalog = catalog))$
#   prep()$
#   bake()
# frec$tidy()
# # ----------------------------------------------------------------------

