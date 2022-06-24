test_that("predict_terms works", {

  library(hydrorecipes)
  library(earthtide)

  data(transducer, package = "hydrorecipes")
  transducer$datetime_num <- as.numeric(transducer$datetime)
  transducer <- transducer[1:5000,]

  rec_toll_rasmussen <- recipe(wl~baro+et+datetime_num, transducer) |>
    step_lead_lag(baro, lag = log_lags(100, 86400 * 2 / 120)) |>
    step_lead_lag(et, lag = seq(0, 180, 20)) |>
    step_ns(datetime_num, deg_free = 3) |>
    prep()
  input_toll_rasmussen <- rec_toll_rasmussen |> bake(new_data = NULL)
  fit_toll_rasmussen <- lm(wl~., input_toll_rasmussen)
  pred <- predict_terms(fit = fit_toll_rasmussen,
                        rec = rec_toll_rasmussen,
                        data = input_toll_rasmussen)

  rec_toll_rasmussen2 <- recipe(wl~baro+et+datetime_num, transducer) |>
    step_lead_lag(baro, lag = log_lags(100, 86400 * 2 / 120), prefix = 'a') |>
    step_lead_lag(et, lag = seq(0, 180, 20), prefix = 'b_') |>
    step_ns(datetime_num, deg_free = 3) |>
    prep()
  input_toll_rasmussen2 <- rec_toll_rasmussen2 |> bake(new_data = NULL)
  fit_toll_rasmussen2 <- lm(wl~., input_toll_rasmussen2)
  pred2 <- predict_terms(fit = fit_toll_rasmussen2,
                        rec = rec_toll_rasmussen2,
                        data = input_toll_rasmussen2)

  expect_equal(pred, pred2)

  expect_equal(nrow(pred), nrow(transducer))
  expect_equal(ncol(pred), 5)
  expect_equal(names(pred), c('lead_lag_baro', 'lead_lag_et', 'ns_datetime_num', 'intercept', 'predicted'))
  expect_equal(sum(is.na(pred)), 86400 * 2 * 2 / 120 + 180)


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
  pred <- predict_terms(fit = fit_rasmussen_mote,
                        rec = rec_rasmussen_mote,
                        data = input_rasmussen_mote)
  expect_equal(nrow(pred), nrow(transducer))
  expect_equal(ncol(pred), 5)
  expect_equal(names(pred), c('lead_lag_baro', 'harmonic_datetime_num', 'ns_datetime_num', 'intercept', 'predicted'))
  expect_equal(sum(is.na(pred)), 86400 * 2 * 2 / 120)



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
  pred <- predict_terms(fit = fit_dl,
                        rec = rec_dl,
                        data = input_dl)
  expect_equal(nrow(pred), nrow(transducer))
  expect_equal(ncol(pred), 5)
  expect_equal(names(pred), c('distributed_lag_baro', 'earthtide_datetime_num', 'ns_datetime_num', 'intercept', 'predicted'))
  expect_equal(sum(is.na(pred)), 86400 * 2 * 2 / 120)


  library(glmnet)
  xy <- na.omit(input_dl)
  x <- as.matrix(input_dl[, -1])
  y <- input_dl[['wl']]
  fit_cv <- cv.glmnet(x, y)
  pred <- predict_terms(fit_cv, rec_dl, xy)
  expect_equal(nrow(pred), nrow(xy))
  expect_equal(ncol(pred), 5)
  expect_equal(names(pred), c('distributed_lag_baro', 'earthtide_datetime_num', 'ns_datetime_num', 'intercept', 'predicted'))
  expect_equal((nrow(transducer) - nrow(pred))*2, 86400 * 2 * 2 / 120)


})
