formula <- as.formula(y~x)
rows <- 20

dat <- data.frame(x = rnorm(rows),
                  y = as.numeric(1:rows),
                  z = rnorm(rows))

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepIntercept$new())$
  plate("tbl")[,c(3,1,2)]
frec[[1]] <- as.integer(frec[[1]])

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_intercept() |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

tinytest::expect_equivalent(frec, rec, info = "frecipes matches recipes")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
frec2 = recipe(formula = formula, data = dat)|>
  step_intercept() |>
  plate("df")

frec1 = Recipe$new(formula = formula, data = dat)$
  add_step(StepIntercept$new())$
  plate("df")

tinytest::expect_equivalent(frec1, frec2,
                            info = "R6 and frecipes api are equivalent")
