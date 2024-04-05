data("kennel_2020")
kennel_2020[, datetime := as.numeric(datetime)]
kennel_2020[, wl2 := wl*0.8]
formula <- as.formula(wl + wl2~.)
n_knots <- 12
deg_free <- 27
max_lag <- 1 + 720

frec = Recipe$new(formula = formula, data = unclass(kennel_2020))$
  add_step(StepDistributedLag$new(baro,
                                  knots = frecipes:::log_lags_arma(n_knots, max_lag)))$
  add_step(StepSplineB$new(datetime, df = deg_free, intercept = FALSE))$
  add_step(StepIntercept$new())$
  add_step(StepDropColumns$new(baro))$
  add_step(StepDropColumns$new(datetime))$
  add_step(StepOlsPredict$new(formula))$
  prep()$
  bake()


