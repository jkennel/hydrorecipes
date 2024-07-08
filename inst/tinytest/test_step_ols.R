data("kennel_2020")
kennel_2020[, datetime := as.numeric(datetime)]
kennel_2020[, wl2 := wl*0.8]
formula <- as.formula(wl + wl2~.)
n_knots <- 12
deg_free <- 27
max_lag <- 1 + 720

formula <- as.formula(wl+wl2~.)
formula2 <- as.formula(wl+wl2~spline_b_datetime_25 )
formula3 <- as.formula(wl+wl2~spline_b_datetime_25 + spline_b_datetime_26)
hrec = hydrorecipes:::Recipe$new(formula = formula, data = unclass(kennel_2020))$
  add_step(hydrorecipes:::StepDistributedLag$new(baro,
                                  knots = hydrorecipes:::log_lags_arma(n_knots, max_lag)))$
  add_step(hydrorecipes:::StepDistributedLag$new(et,
                                  knots = hydrorecipes:::log_lags_arma(n_knots, max_lag)))$
  add_step(hydrorecipes:::StepSplineB$new(datetime, df = deg_free, intercept = FALSE))$
  add_step(hydrorecipes:::StepIntercept$new())$
  add_step(hydrorecipes:::StepDropColumns$new(baro))$
  # add_step(hydrorecipes:::StepDropColumns$new(et))$
  add_step(hydrorecipes:::StepDropColumns$new(datetime))$
  add_step(hydrorecipes:::StepOls$new(formula))$
  prep()$
  bake()

hrec$get_response_data(type = 'dt')
