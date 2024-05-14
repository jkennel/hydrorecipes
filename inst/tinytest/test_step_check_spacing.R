formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = rnorm(rows))


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# hydrorecipes version
frec_irr = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepCheckSpacing$new(x))$
  prep()$
  bake()$
  get_step_data("check")
expect_equivalent(frec_irr[[1]], FALSE,
                            info = "irregular spacing")


frec_reg = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepCheckSpacing$new(y))$
  prep()$
  bake()$
  get_step_data("check")
expect_equivalent(frec_reg[[1]], TRUE,
                            info = "regular spacing")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

frec1 = recipe(formula = formula, data = dat) |>
  step_check_spacing(x) |>
  prep() |>
  bake()

frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepCheckSpacing$new(x))$prep()$bake()


expect_equivalent(frec1$get_step_data("check"),
                            frec2$get_step_data("check"),
                            info = "R6 and hydrorecipes api are equivalent")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
