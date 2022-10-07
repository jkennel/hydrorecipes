test_that("step_find_interval works", {
  data(wipp30)

  breaks <- c(9000, 12000)

  r_fi_one_hot <- recipe(wl~., data = wipp30) |>
    step_find_interval(time,
                       vec = breaks,
                       keep_original_cols = TRUE,
                       encoding = 'one_hot') |>
    prep() |>
    bake(new_data = NULL)

  r_fi_dummy <- recipe(wl~., data = wipp30) |>
    step_find_interval(time,
                       vec = breaks,
                       keep_original_cols = TRUE,
                       encoding = 'dummy') |>
    prep() |>
    bake(new_data = NULL)

  r_fact <- recipe(wl~., data = wipp30) |>
    step_find_interval(time, vec = breaks, encoding = 'factor') |>
    prep() |>
    bake(new_data = NULL)

  r_int <- recipe(wl~., data = wipp30) |>
    step_find_interval(time, vec = breaks, encoding = 'integer') |>
    prep() |>
    bake(new_data = NULL)

  wipp30$intervals <- factor(findInterval(wipp30$time, vec = breaks))

  r_dummy <- recipe(wl~., data = wipp30) |>
    step_dummy(intervals) |>
    prep() |>
    bake(new_data = NULL)

  r_one_hot <- recipe(wl~., data = wipp30) |>
    step_dummy(intervals, one_hot = TRUE) |>
    prep() |>
    bake(new_data = NULL)


  expect_equal(r_int$find_interval_time,
               findInterval(wipp30$time, vec = breaks),
               ignore_attr = TRUE)
  expect_equal(r_fact$find_interval_time, wipp30$intervals)
  expect_equal(r_dummy, r_fi_dummy, ignore_attr = TRUE)
  expect_equal(r_one_hot, r_fi_one_hot, ignore_attr = TRUE)

  expect_error(recipe(wl~., data = wipp30) |>
                 step_find_interval(time))

})
