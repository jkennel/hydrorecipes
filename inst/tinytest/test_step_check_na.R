formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = rnorm(rows))


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# frecipes version
frec_false = Recipe$new(formula = formula, data = dat)$
  add_step(StepCheckNA$new(x))$
  prep()$
  bake()$
  get_step_data("check")
tinytest::expect_equivalent(frec_false[[1]], FALSE,
                            info = "No NAs present")


dat[10,1:3] <- NA_real_
frec_true = Recipe$new(formula = formula, data = dat)$
  add_step(StepCheckNA$new(y))$
  prep()$
  bake()$
  get_step_data("check")
tinytest::expect_equivalent(frec_true, TRUE,
                            info = "NAs present")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

frec1 = recipe(formula = formula, data = dat) |>
  step_check_na(x) |>
  prep() |>
  bake()

frec2 = Recipe$new(formula = formula, data = dat)$
  add_step(StepCheckNA$new(x))$bake()


tinytest::expect_equivalent(frec1$get_step_data("check"), frec2$get_step_data("check"),
                            info = "R6 and frecipes api are equivalent")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
