test_that("step_distributed_lag works", {

  data(wipp30)
  max_time_lag <- 1000
  n_lags <- 10
  rec <- recipe(wl~baro, wipp30)


  expect_equal(rec |>
                 step_distributed_lag(baro,
                                      knots = c(0, 6000)) |>
                 prep() |>
                 bake(new_data = NULL),
               rec |>
                 step_distributed_lag(baro,
                                      knots = c(0,6000),
                                      basis_mat = basis_lag(c(0,6000))) |>
                 prep() |>
                 bake(new_data = NULL)
  )
  m <- as.matrix(basis_lag(c(0,6000)))
  attr(m, 'knots') <- NULL
  expect_equal(rec |>
                 step_distributed_lag(baro,
                                      knots = c(0, 6000)) |>
                 prep() |>
                 bake(new_data = NULL),
               rec |>
                 step_distributed_lag(baro,
                                      knots = c(0,6000),
                                      basis_mat = m) |>
                 prep() |>
                 bake(new_data = NULL)
  )


  expect_error(rec |>
                 step_distributed_lag(baro,
                                      knots = c(0,5000),
                                      basis_mat = m) |>
                 prep() |>
                 bake(new_data = NULL)
  )

  expect_s3_class(rec |>
                    step_distributed_lag(baro, knots = log_lags(10, max_time_lag)) |>
                    prep() |>
                    bake(new_data = NULL), 'data.frame')

  expect_s3_class(rec |>
                    step_distributed_lag(baro, knots = c(0, 6000)) |>
                    prep() |>
                    bake(new_data = NULL), 'data.frame')


  expect_equal(nrow(na.omit(rec |>
                              step_distributed_lag(baro, knots = log_lags(n_lags, max_time_lag)) |>
                              prep() |>
                              bake(new_data = NULL))), nrow(wipp30) - max_time_lag)

  expect_equal(ncol(rec |>
                      step_distributed_lag(baro, knots = log_lags(n_lags, max_time_lag)) |>
                      prep() |>
                      bake(new_data = NULL)), n_lags + 1)

  expect_equal(ncol(rec |>
                      step_distributed_lag(baro, knots = log_lags(n_lags, max_time_lag), keep_original_cols = TRUE) |>
                      prep() |>
                      bake(new_data = NULL)), n_lags + 2)



  # expect_warning(rec |>
  #                  step_distributed_lag(baro, knots = c(1,1, 2)) |>
  #                  prep() |>
  #                  bake(new_data = NULL))
  # expect_warning(rec |>
  #                  step_distributed_lag(baro, knots = c(1, 1, 2)) |>
  #                  prep() |>
  #                  bake(new_data = NULL))

  expect_error(rec |>
                 step_distributed_lag(baro, knots = 1) |>
                 prep() |>
                 bake(new_data = NULL))
  expect_error(rec |>
                 step_distributed_lag(baro, knots = -10:(-1)) |>
                 prep() |>
                 bake(new_data = NULL))
  expect_error(rec |>
                 step_distributed_lag(baro, knots = c(0, nrow(wipp30) + 2)) |>
                 prep() |>
                 bake(new_data = NULL))
  expect_error(rec |>
                 step_distributed_lag(baro, knots = c(1, 1)) |>
                 prep() |>
                 bake(new_data = NULL))

  rec_prep <- rec |>
    step_distributed_lag(baro, knots = log_lags(10, max_time_lag)) |>
    prep()
  expect_equal(tidy(rec_prep, 1), tidy2(rec_prep, 1))

  expect_output(print(rec_prep, 1))

  rec2 <- rec |>
    step_distributed_lag(baro, knots = log_lags(10, max_time_lag))
  expect_equal(tidy2(rec_prep, 1)$key, tidy2(rec2, 1)$key)


  # test different spline models
  expect_equal(r1 <- rec |>
                 step_distributed_lag(baro,
                                      spline_fun = splines::ns,
                                      options = list(intercept = TRUE),
                                      knots = c(0, 6000)) |>
                 prep() |>
                 bake(new_data = NULL),
               r2 <- rec |>
                 step_distributed_lag(baro,
                                      knots = c(0, 6000)) |>
                 prep() |>
                 bake(new_data = NULL)
  )


})
