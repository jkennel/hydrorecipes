test_that("tidy2 works", {

  data(wipp30)
  max_time_lag <- 1000
  n_lags <- 10
  trained <- recipe(wl~., wipp30) |>
                         step_distributed_lag(baro, knots = log_lags(10, max_time_lag)) |>
                         prep()
  expect_s3_class(tidy2(trained, 1), 'data.frame')

  expect_s3_class(tidy2(trained), 'data.frame')

  expect_error(tidy2(trained, number = NULL))
  expect_error(tidy2(trained, number = 100))
  expect_error(tidy2(trained, number = 1, id = "id"))
  expect_error(tidy2(trained, id = "id"))

  trained <- recipe(wl~., wipp30) |>
    step_ns(baro, deg_free = 5) |>
    prep()
  expect_s3_class(tidy2(trained, 1), 'data.frame')

  trained <- recipe(wl~., wipp30) |>
    step_intercept() |>
    prep()

  expect_s3_class(tidy(trained, 1), 'data.frame')
  expect_s3_class(tidy2(trained, 1), 'data.frame')
  expect_equal(tidy2(trained, 1), tidy(trained, 1))



  trained <- recipe(wl~., wipp30) |>
    step_ns()

  expect_s3_class(tidy2(trained), 'data.frame')

})

