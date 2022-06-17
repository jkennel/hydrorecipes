test_that("step_lead_lag works", {
  lags <- -5:5
  n_lags <- length(lags)

  df <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, keep_original_cols = TRUE) |>
    prep() |>
    bake(new_data = NULL)

  expect_s3_class(df, 'data.frame')
  expect_equal(nrow(na.omit(df)), nrow(wipp30) - n_lags + 1)
  expect_equal(ncol(df), ncol(wipp30) + n_lags)


  df <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL)

  expect_equal(ncol(df), ncol(wipp30) + n_lags - 1)

  expect_warning(recipe(wl~., wipp30) |>
                      step_lead_lag(baro, lag = c(1,1,2)) |>
                      prep() |>
                      bake(new_data = NULL))


  df <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, n_subset = 2, keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL)

  expect_equal(nrow(df), ceiling(nrow(wipp30) / 2))

  expect_error(recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, n_subset = 0, keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL))


  # n_subset
  a <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, n_subset = 2, keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL)
  b <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, n_subset = 1, keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL)
  expect_equal(a, b[seq(1, nrow(wipp30), 2),])


  a <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, n_subset = 3, keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL)
  b <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, n_subset = 1, keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL)
  expect_equal(a, b[seq(1, nrow(wipp30), 3),])

  expect_error(recipe(wl~., wipp30) |>
                 step_lead_lag(baro, lag = 100000, keep_original_cols = FALSE) |>
                 prep() |>
                 bake(new_data = NULL))


  # n_shift

  a <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, n_subset = 2, n_shift = 1, keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL)
  b <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, n_subset = 1, n_shift = 0, keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL)

  expect_equal(a, b[seq(2, nrow(wipp30), 2),])

  a <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, n_subset = 3, n_shift = 2, keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL)
  b <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, n_subset = 1, n_shift = 0, keep_original_cols = FALSE) |>
    prep() |>
    bake(new_data = NULL)

  expect_equal(a, b[seq(3, nrow(wipp30), 3),])

  expect_error(recipe(wl~., wipp30) |>
                   step_lead_lag(baro, lag = lags, n_subset = 1, n_shift = 1, keep_original_cols = FALSE) |>
                   prep() |>
                   bake(new_data = NULL))
  expect_error(recipe(wl~., wipp30) |>
                 step_lead_lag(baro, lag = lags, n_subset = 2, n_shift = 2, keep_original_cols = FALSE) |>
                 prep() |>
                 bake(new_data = NULL))

  rec <- recipe(wl~., wipp30) |>
    step_lead_lag(baro, lag = lags, n_subset = 3, n_shift = 2, keep_original_cols = FALSE) |>
    prep()

  expect_equal(tidy(rec, 1), tidy2(rec, 1))


})
