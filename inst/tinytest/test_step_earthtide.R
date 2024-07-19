data(kennel_2020)

# synthetic earthtide
latitude     <- 34.23411                           # latitude
longitude    <- -118.678                           # longitude
elevation    <- 500                                # elevation
cutoff       <- 1e-5                               # cutoff
catalog      <- 'ksm04'                            # hartmann wenzel catalog
astro_update <- 300                                # how often to update astro parameters
method       <- 'volume_strain'                    # which potential to calculate

wave_groups_dl <- earthtide::eterna_wavegroups
wave_groups_dl <- na.omit(wave_groups_dl[wave_groups_dl$time == '1 month', ])
wave_groups_dl <- wave_groups_dl[wave_groups_dl$start > 0.5,]
wave_groups_dl <- wave_groups_dl[, c("start", "end")]
ngr <- nrow(wave_groups_dl)


frec = hydrorecipes:::Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
  add_step(hydrorecipes:::StepEarthtide$new(datetime,
                                        wave_groups = wave_groups_dl,
                                        latitude = latitude,
                                        longitude = longitude,
                                        elevation = elevation,
                                        cutoff = cutoff,
                                        catalog = catalog))$
  plate()
et <- earthtide::calc_earthtide(kennel_2020$datetime,
                                wave_groups = wave_groups_dl,
                                latitude = latitude,
                                longitude = longitude,
                                elevation = elevation,
                                cutoff = cutoff,
                                catalog = catalog)

expect_equivalent(frec$earthtide, et$gravity)


frec = hydrorecipes:::Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
  add_step(hydrorecipes:::StepEarthtide$new(datetime,
                                        do_predict = FALSE,
                                        wave_groups = wave_groups_dl,
                                        latitude = latitude,
                                        longitude = longitude,
                                        elevation = elevation,
                                        cutoff = cutoff,
                                        catalog = catalog))$
  plate()
et <- earthtide::calc_earthtide(kennel_2020$datetime,
                                do_predict = FALSE,
                                wave_groups = wave_groups_dl,
                                latitude = latitude,
                                longitude = longitude,
                                elevation = elevation,
                                cutoff = cutoff,
                                catalog = catalog)

expect_equivalent(frec[,-c(2L, 3L)], et)

# bench::mark(
# frec1 = hydrorecipes:::Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
#   add_step(hydrorecipes:::StepEarthtide$new(datetime,
#                                             do_predict = TRUE,
#                                             wave_groups = wave_groups_dl,
#                                             latitude = latitude,
#                                             longitude = longitude,
#                                             elevation = elevation,
#                                             cutoff = cutoff,
#                                             catalog = catalog,
#                                             interp_factor = 10L,
#                                             n_thread = 10L))$
#   plate(),
#
# frec2 = hydrorecipes:::Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
#   add_step(hydrorecipes:::StepEarthtide$new(datetime,
#                                             do_predict = TRUE,
#                                             wave_groups = wave_groups_dl,
#                                             latitude = latitude,
#                                             longitude = longitude,
#                                             elevation = elevation,
#                                             cutoff = cutoff,
#                                             catalog = catalog,
#                                             interp_factor = 1L))$
#   plate(), check = FALSE
# )
#
# plot(earthtide~datetime, frec1, type = "l")
# plot(earthtide~datetime, frec2, type = "l", col = "red")
