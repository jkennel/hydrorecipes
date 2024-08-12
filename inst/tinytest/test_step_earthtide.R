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


hrec = hydrorecipes:::Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
  add_step(hydrorecipes:::StepEarthtide$new(datetime,
                                        wave_groups = wave_groups_dl,
                                        latitude = latitude,
                                        longitude = longitude,
                                        elevation = elevation,
                                        cutoff = cutoff,
                                        catalog = catalog,
                                        astro_update = 60))$
  plate()
et <- earthtide::calc_earthtide(kennel_2020$datetime,
                                wave_groups = wave_groups_dl,
                                latitude = latitude,
                                longitude = longitude,
                                elevation = elevation,
                                cutoff = cutoff,
                                catalog = catalog,
                                astro_update = 60)

expect_equivalent(hrec$earthtide, et$gravity)


hrec = hydrorecipes:::Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
  add_step(hydrorecipes:::StepEarthtide$new(datetime,
                                        do_predict = FALSE,
                                        wave_groups = wave_groups_dl,
                                        latitude = latitude,
                                        longitude = longitude,
                                        elevation = elevation,
                                        cutoff = cutoff,
                                        catalog = catalog,
                                        astro_update = 60))$
  plate()
et <- earthtide::calc_earthtide(kennel_2020$datetime,
                                do_predict = FALSE,
                                wave_groups = wave_groups_dl,
                                latitude = latitude,
                                longitude = longitude,
                                elevation = elevation,
                                cutoff = cutoff,
                                catalog = catalog,
                                astro_update = 60)

hrec2 = recipe(formula = wl~baro+datetime, data = kennel_2020) |>
  step_earthtide(datetime,
                 do_predict = FALSE,
                 wave_groups = wave_groups_dl,
                 latitude = latitude,
                 longitude = longitude,
                 elevation = elevation,
                 cutoff = cutoff,
                 catalog = catalog,
                 astro_update = 60) |>
  plate()

expect_equivalent(hrec, hrec2)
expect_equivalent(hrec[,-c(2L, 3L)], et)

hrec1 = hydrorecipes:::Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
  add_step(hydrorecipes:::StepEarthtide$new(datetime,
                                            do_predict = TRUE,
                                            wave_groups = wave_groups_dl,
                                            latitude = latitude,
                                            longitude = longitude,
                                            elevation = elevation,
                                            cutoff = cutoff,
                                            catalog = catalog,
                                            interp_factor = 5L,
                                            n_thread = 10L))$
  plate()

hrec2 = hydrorecipes:::Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
  add_step(hydrorecipes:::StepEarthtide$new(datetime,
                                            do_predict = TRUE,
                                            wave_groups = wave_groups_dl,
                                            latitude = latitude,
                                            longitude = longitude,
                                            elevation = elevation,
                                            cutoff = cutoff,
                                            catalog = catalog,
                                            interp_factor = 1L,
                                            n_thread = 10L))$
  plate()
expect_equivalent(hrec1, hrec2)


# plot(earthtide~datetime, hrec1, type = "l")
# plot(earthtide~datetime, hrec2, type = "l", col = "red")
