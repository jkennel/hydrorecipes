data("kennel_2020")
kennel_2020[, datetime := as.numeric(datetime)]
kennel_2020[, wl2 := wl * 0.8]
kennel_2020[, wl3 := wl * 0.6]
formula <- as.formula(wl + wl2 + wl3~.)
n_knots <- 12
deg_free <- 27
lag_max <- 1 + 720

formula <- as.formula(wl+wl2 + wl3~.)
formula2 <- as.formula(wl+wl2~spline_b_datetime_25 )
formula3 <- as.formula(wl+wl2~spline_b_datetime_25 + spline_b_datetime_26)
formula4 <- as.formula(wl+wl2+wl3~spline_b_datetime_25 + spline_b_datetime_26)
hrec = hydrorecipes:::Recipe$new(formula = formula, data = unclass(kennel_2020))$
  add_step(hydrorecipes:::StepDistributedLag$new(baro,
                                  knots = hydrorecipes:::log_lags(n_knots, lag_max)))$
  add_step(hydrorecipes:::StepDistributedLag$new(et,
                                  knots = hydrorecipes:::log_lags(n_knots, lag_max)))$
  add_step(hydrorecipes:::StepSplineB$new(datetime, df = deg_free, intercept = FALSE))$
  add_step(hydrorecipes:::StepIntercept$new())$
  add_step(hydrorecipes:::StepDropColumns$new(baro))$
  add_step(hydrorecipes:::StepHarmonic$new(datetime, frequency = c(1, 2, 3), cycle_size = 86400))$
  # add_step(hydrorecipes:::StepEarthtide$new(datetime, do_predict = FALSE, astro_update = 30))$
  # add_step(hydrorecipes:::StepDropColumns$new(et))$
  add_step(hydrorecipes:::StepDropColumns$new(datetime))$
  add_step(hydrorecipes:::StepOls$new(formula))$
  add_step(hydrorecipes:::StepOls$new(formula2))$
  add_step(hydrorecipes:::StepOls$new(formula3))$
  add_step(hydrorecipes:::StepOls$new(formula4))$
  prep()$
  bake()

hrec$get_response_data(type = 'dt')[grep("harmonic", step_id)]
hrec$get_response_data(type = 'dt')[grep("earthtide", step_id)]
hrec$get_response_data(type = 'dt')[grep("spline", step_id)]





hrec = hydrorecipes:::Recipe$new(formula = formula, data = unclass(kennel_2020))$
  add_step(hydrorecipes:::StepDistributedLag$new(baro,
                                                 knots = hydrorecipes:::log_lags(n_knots, lag_max)))$
  add_step(hydrorecipes:::StepDistributedLag$new(et,
                                                 knots = hydrorecipes:::log_lags(n_knots, lag_max)))$
  add_step(hydrorecipes:::StepSplineB$new(datetime, df = deg_free, intercept = FALSE))$
  add_step(hydrorecipes:::StepIntercept$new())$
  add_step(hydrorecipes:::StepDropColumns$new(baro))$
  add_step(hydrorecipes:::StepHarmonic$new(datetime, frequency = c(1, 2, 3), cycle_size = 86400))$
  add_step(hydrorecipes:::StepEarthtide$new(datetime, do_predict = FALSE, astro_update = 30))$
  # add_step(hydrorecipes:::StepDropColumns$new(et))$
  add_step(hydrorecipes:::StepDropColumns$new(datetime))$
  add_step(hydrorecipes:::StepOls$new(formula))$
  prep()$
  bake()

expect_equivalent(class(hrec$get_response_data(type = 'dt')),
                  c("data.table", "data.frame"))
expect_equivalent(nrow(hrec$get_response_data(type = 'dt')),
                  5850L)
expect_equivalent(class(hrec$get_response_data(type = 'df')),
                  "data.frame")
