test_that("multiplication works", {
  n <- 100
  library(earthtide)
  data(eterna_wavegroups)

  wg <- na.omit(eterna_wavegroups[eterna_wavegroups$time == '1 month',])

  d <- data.frame(fct = as.factor(sample(5, n, replace = TRUE)),
                  num = rnorm(n),
                  num2 = rnorm(n),
                  num3 = rnorm(n),
                  cnt = 1:n,
                  dd = as.POSIXct(1:n, origin = '2021-01-01'),
                  out = rnorm(n)
  )

  r <-  recipe(out~., d) |>
    step_harmonic(cnt, frequency = c(20), cycle_size = 10, keep_original_cols = TRUE) |>
    step_earthtide(dd,
                   latitude = 34,
                   longitude = -118.5,
                   wave_groups = wg,
                   do_predict = FALSE) |>
    step_ns(num, deg_free = 2, ) |>
    step_dummy(fct, one_hot = TRUE) |>
    step_lag(num3, lag = 0:2) |>
    step_lead_lag(num2, lag = 0:2) |>
    prep()

  rec_names <- unlist(recipe_term_names(r))

  inter <- intersect(unlist(recipe_term_names(r)),
            names(r |> bake(new_data = NULL)))

  expect_true(length(rec_names)==length(inter))

  r <-  recipe(out~., d) |>
    step_harmonic(cnt, frequency = c(20), cycle_size = 10, keep_original_cols = TRUE) |>
    step_earthtide(dd,
                   latitude = 34,
                   longitude = -118.5,
                   wave_groups = wg,
                   do_predict = FALSE) |>
    step_ns(num, deg_free = 2, ) |>
    step_dummy(fct) |>
    step_lag(num3, lag = 0:2) |>
    step_lead_lag(num2, lag = 0:2) |>
    prep()

  rec_names <- unlist(recipe_term_names(r))

  inter <- intersect(unlist(recipe_term_names(r)),
                     names(r |> bake(new_data = NULL)))

  expect_true(length(rec_names)==length(inter))


})
