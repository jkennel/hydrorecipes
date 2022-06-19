test_that("step_distributed_lag works", {

  data(wipp30)
  max_time_lag <- 1000
  n_lags <- 10
  expect_s3_class(recipe(wl~., wipp30) |>
                         step_distributed_lag(baro, knots = log_lags(10, max_time_lag)) |>
                         prep() |>
                         bake(new_data = NULL), 'data.frame')

  expect_s3_class(recipe(wl~., wipp30) |>
                         step_distributed_lag(baro, knots = c(0, 6000)) |>
                         prep() |>
                         bake(new_data = NULL), 'data.frame')


  expect_equal(nrow(na.omit(recipe(wl~., wipp30) |>
                 step_distributed_lag(baro, knots = log_lags(n_lags, max_time_lag)) |>
                 prep() |>
                 bake(new_data = NULL))), nrow(wipp30) - max_time_lag)

  expect_equal(ncol(recipe(wl~., wipp30) |>
                              step_distributed_lag(baro, knots = log_lags(n_lags, max_time_lag)) |>
                              prep() |>
                              bake(new_data = NULL)), ncol(wipp30) + n_lags - 1)

  expect_equal(ncol(recipe(wl~., wipp30) |>
                      step_distributed_lag(baro, knots = log_lags(n_lags, max_time_lag), keep_original_cols = TRUE) |>
                      prep() |>
                      bake(new_data = NULL)), ncol(wipp30) + n_lags)



  expect_warning(recipe(wl~., wipp30) |>
                              step_distributed_lag(baro, knots = c(1,1, 2)) |>
                              prep() |>
                              bake(new_data = NULL))
  expect_warning(recipe(wl~., wipp30) |>
                 step_distributed_lag(baro, knots = c(1, 1, 2)) |>
                 prep() |>
                 bake(new_data = NULL))

  expect_error(recipe(wl~., wipp30) |>
                 step_distributed_lag(baro, knots = 1) |>
                 prep() |>
                 bake(new_data = NULL))
  expect_error(recipe(wl~., wipp30) |>
                 step_distributed_lag(baro, knots = -10:(-1)) |>
                 prep() |>
                 bake(new_data = NULL))
  expect_error(recipe(wl~., wipp30) |>
                 step_distributed_lag(baro, knots = c(0, nrow(wipp30) + 2)) |>
                 prep() |>
                 bake(new_data = NULL))
  expect_error(recipe(wl~., wipp30) |>
                   step_distributed_lag(baro, knots = c(1, 1)) |>
                   prep() |>
                   bake(new_data = NULL))

  rec <- recipe(wl~., wipp30) |>
    step_distributed_lag(baro, knots = log_lags(10, max_time_lag)) |>
    prep()
  expect_equal(tidy(rec, 1), tidy2(rec, 1))

  expect_output(print(rec, 1))

  rec2 <- recipe(wl~., wipp30) |>
    step_distributed_lag(baro, knots = log_lags(10, max_time_lag))
  expect_equal(tidy2(rec, 1)$key, tidy2(rec2, 1)$key)

})
