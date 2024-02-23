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
  checks
tinytest::expect_equivalent(frec_irr[[1]], FALSE,
                            info = "irregular spacing")


dat[10,1] <- NA_real_
frec_true = Recipe$new(formula = formula, data = dat)$
  add_step(StepCheckNA$new(y))$
  prep()$
  bake()$
  checks
tinytest::expect_equivalent(frec_reg[[1]], TRUE,
                            info = "regular spacing")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

frec1 = recipe(formula = formula, data = dat) |>
  step_check_na(x) |>
  bake()

frec2 = Recipe$new(formula = formula, data = dat)$
  add_step(StepCheckNA$new(x))$bake()$checks


tinytest::expect_equivalent(frec1$checks, frec2,
                            info = "R6 and frecipes api are equivalent")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
