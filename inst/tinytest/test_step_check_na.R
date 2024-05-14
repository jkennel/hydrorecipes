formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = rnorm(rows))


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# hydrorecipes version
frec_false = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepCheckNA$new(x))$
  prep()$
  bake()$
  get_step_data("check")
expect_equivalent(unlist(frec_false[[1]][["x"]]), FALSE,
                            info = "No NAs present")


dat[10,1:3] <- NA_real_
frec_true = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepCheckNA$new(y))$
  prep()$
  bake()$
  get_step_data("check")
expect_equivalent(unlist(frec_true[[1]][["y"]]), TRUE,
                            info = "NAs present")

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

frec1 = recipe(formula = formula, data = dat) |>
  step_check_na(x) |>
  prep() |>
  bake()

frec2 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepCheckNA$new(x))$bake()


expect_equivalent(frec1$get_step_data("check"), frec2$get_step_data("check"),
                            info = "R6 and hydrorecipes api are equivalent")
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
