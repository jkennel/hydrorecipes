test_that("response works", {

  library(hydrorecipes)
  library(earthtide)

  data(transducer, package = "hydrorecipes")
  transducer <- transducer[1:5000,]
  transducer$datetime_num <- as.numeric(transducer$datetime)

  rec_toll_rasmussen <- recipe(wl~baro+et+datetime_num, transducer) |>
    step_lead_lag(baro, lag = log_lags(100, 86400 * 2 / 120)) |>
    step_lead_lag(et, lag = seq(0, 180, 20)) |>
    step_ns(datetime_num, deg_free = 3) |>
    prep()
  input_toll_rasmussen <- rec_toll_rasmussen |> bake(new_data = NULL)
  fit_toll_rasmussen <- lm(wl~., input_toll_rasmussen)
  resp <- response(fit_toll_rasmussen, rec_toll_rasmussen)
  expect_equal(length(unique(resp$term)), 2)
  expect_equal(nrow(resp), 200+20)
  expect_equal(unique(resp[resp$term == 'baro',]$x), log_lags(100, 86400 * 2 / 120))
  expect_equal(unique(resp[resp$term == 'et',]$x), seq(0, 180, 20))


  tidal_freqs <- c(0.89324406, 0.92953571, 0.96644626, 1.00273791, 1.03902956,
                   1.07594011, 1.86454723, 1.89598197, 1.93227362, 1.96856526,
                   2.00000000, 2.89841042)
  rec_rasmussen_mote <- recipe(wl~baro+datetime_num, transducer) |>
    step_lead_lag(baro, lag = log_lags(100, 86400 * 2 / 120)) |>
    step_harmonic(datetime_num,
                  frequency = tidal_freqs,
                  cycle_size = 86400,
                  keep_original_cols = TRUE) |>
    step_ns(datetime_num, deg_free = 3) |>
    prep()
  input_rasmussen_mote <- rec_rasmussen_mote |> bake(new_data = NULL)
  fit_rasmussen_mote <- lm(wl~., input_rasmussen_mote)
  resp <- response(fit_rasmussen_mote, rec_rasmussen_mote)
  expect_equal(length(unique(resp$term)), 2)
  expect_equal(nrow(resp), 200+length(tidal_freqs) * 4)
  expect_equal(unique(resp[resp$term == 'baro',]$x), log_lags(100, 86400 * 2 / 120))
  expect_equal(unique(resp[resp$term == 'datetime_num',]$x), tidal_freqs)



  wave_groups <- earthtide::eterna_wavegroups
  wave_groups <- na.omit(wave_groups[wave_groups$time == '1 month', ])
  wave_groups <- wave_groups[wave_groups$start > 0.5, ]
  latitude  <- 34.0
  longitude <- -118.5
  rec_dl <- recipe(wl~baro+datetime_num, transducer) |>
    step_distributed_lag(baro, knots = log_lags(15, 86400 * 2 / 120)) |>
    step_earthtide(datetime_num,
                   latitude = latitude,
                   longitude = longitude,
                   astro_update = 1,
                   wave_groups = wave_groups) |>
    step_ns(datetime_num, deg_free = 3) |>
    prep()

  input_dl <- rec_dl |> bake(new_data = NULL)
  fit_dl <- lm(wl~., input_dl)
  resp <- response(fit_dl, rec_dl)
  expect_equal(length(unique(resp$term)), 2)
  expect_equal(nrow(resp), 86400 * 2 / 120 * 2 + 2 + nrow(wave_groups) * 4)
  expect_equal(unique(resp[resp$term == 'baro',]$x), 0:(86400*2/120))
  expect_equal(unique(resp[resp$term == 'datetime_num',]$x), tidal_freqs)

  library(glmnet)
  xy <- na.omit(input_dl)
  x <- as.matrix(xy[, -1])
  y <- xy[['wl']]
  fit_cv <- cv.glmnet(x, y, family = 'gaussian', alpha = 0.1)
  resp <- response(fit_cv, rec_dl)
  expect_equal(length(unique(resp$term)), 2)
  expect_equal(nrow(resp), 86400 * 2 / 120 * 2 + 2 + nrow(wave_groups) * 4)
  expect_equal(unique(resp[resp$term == 'baro',]$x), 0:(86400*2/120))
  expect_equal(unique(resp[resp$term == 'datetime_num',]$x), tidal_freqs)


  expect_output(tmp <- response(fit_dl, rec_dl, verbose = TRUE))


})
