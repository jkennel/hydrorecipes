data(kennel_2020)

# synthetic earthtide
latitude     <- 34.23411                           # latitude
longitude    <- -118.678                           # longitude
elevation    <- 500                                # elevation
cutoff       <- 1e-5                               # cutoff
catalog      <- 'ksm04'                            # hartmann wenzel catalog
astro_update <- 300                                # how often to update astro parameters
method       <- 'volume_strain'                    # which potential to calculate

wave_groups_dl <- as.data.table(earthtide::eterna_wavegroups)
wave_groups_dl <- na.omit(wave_groups_dl[time == '1 month'])
wave_groups_dl <- wave_groups_dl[wave_groups_dl$start > 0.5,]
wave_groups_dl <- wave_groups_dl[, list(start, end)]
ngr <- nrow(wave_groups_dl)


frec = Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
  add_step(StepEarthtide$new(datetime,
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

tinytest::expect_equivalent(frec$earthtide, et$gravity)


frec = Recipe$new(formula = wl~baro+datetime, data = kennel_2020)$
  add_step(StepEarthtide$new(datetime,
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

tinytest::expect_equivalent(frec[,-c(2L, 3L)], et)

