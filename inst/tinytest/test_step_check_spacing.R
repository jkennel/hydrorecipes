formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = rnorm(rows))


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# frecipes version
frec_irr = Recipe$new(formula = formula, data = dat)$
  add_step(StepCheckSpacing$new(x))$
  prep()$
  bake()$
  get_step_data("check")
tinytest::expect_equivalent(frec_irr[[1]], FALSE,
                            info = "irregular spacing")


frec_reg = Recipe$new(formula = formula, data = dat)$
  add_step(StepCheckSpacing$new(y))$
  prep()$
  bake()$
  get_step_data("check")
tinytest::expect_equivalent(frec_reg[[1]], TRUE,
                            info = "regular spacing")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

frec1 = recipe(formula = formula, data = dat) |>
  step_check_spacing(x) |>
  prep() |>
  bake()

frec2 = Recipe$new(formula = formula, data = dat)$
  add_step(StepCheckSpacing$new(x))$prep()$bake()


tinytest::expect_equivalent(frec1$get_step_data("check"),
                            frec2$get_step_data("check"),
                            info = "R6 and frecipes api are equivalent")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
