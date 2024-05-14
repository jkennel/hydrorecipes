set.seed(1)

formula <- as.formula(x~a+b+d+e+f+g)
rows <- 1000

    dat <- data.frame(x = rnorm(rows),
                      a = rnorm(rows),
                      b = rnorm(rows),
                      d = rnorm(rows),
                      e = rnorm(rows),
                      f = rnorm(rows),
                      g = rnorm(rows))
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# hydrorecipes version
frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepPca$new(all_numeric(), n_comp = 6, scale = TRUE, center = TRUE))$
  plate("tbl")
# recipes version
rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_scale(recipes::all_numeric()) |>
  recipes::step_center(recipes::all_numeric()) |>
  recipes::step_pca(recipes::all_numeric(), num_comp = 6, keep_original_cols = FALSE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

expect_equivalent(abs(frec[, -(1:ncol(dat))]), abs(rec))

frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepPca$new(all_numeric(), n_comp = 6, scale = FALSE, center = FALSE))$
  plate("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_pca(recipes::all_numeric(), num_comp = 6, keep_original_cols = FALSE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

expect_equivalent(abs(frec[, -(1:ncol(dat))]), abs(rec))

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# R6 version
frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepPca$new(all_numeric()))$
  plate("df")
# standard version
rec  = recipe(formula = formula, data = dat) |>
  step_pca(all_numeric()) |>
  plate()

expect_equivalent(frec, rec, info = "StepPca with recipes api")

