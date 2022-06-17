test_that("step_earthtide works", {
  library(earthtide)
  data(eterna_wavegroups)

  data(transducer)
  wg <- na.omit(eterna_wavegroups[eterna_wavegroups$time == '1 month',])
  t_sub <- transducer[1:1000,]

  expect_equal(ncol(recipe(wl ~ ., data = t_sub) |>
                      step_earthtide(datetime,
                                     latitude = 34,
                                     longitude = -118.5,
                                     wave_groups = wg,
                                     do_predict = FALSE) |>
                      prep() |>
                      bake(new_data = NULL)), ncol(t_sub) + 2*nrow(wg))

  expect_equal(ncol(recipe(wl ~ ., data = t_sub) |>
                      step_earthtide(datetime,
                                     latitude = 34,
                                     longitude = -118.5,
                                     wave_groups = wg,
                                     do_predict = TRUE) |>
                      prep() |>
                      bake(new_data = NULL)), ncol(t_sub) + 1)

  expect_true(any(names(recipe(wl ~ ., data = t_sub) |>
                          step_earthtide(datetime,
                                         latitude = 34,
                                         longitude = -118.5,
                                         wave_groups = wg,
                                         do_predict = TRUE) |>
                          prep() |>
                          bake(new_data = NULL)) %in% "earthtide_datetime"))

  rec <- recipe(wl ~ ., data = t_sub) |>
    step_earthtide(datetime,
                   latitude = 34,
                   longitude = -118.5,
                   wave_groups = wg,
                   do_predict = TRUE) |>
    prep()

  expect_equal(tidy(rec, 1), tidy2(rec, 1))
  expect_output(print(rec))



})
