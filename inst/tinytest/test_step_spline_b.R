set.seed(1)

formula <- as.formula(x~y+z)
rows <- 1e5

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = cumsum(rnorm(rows)))
ik <- collapse::fquantile(dat$x, probs = seq(0, 1, 0.1))
bk <- ik[c(1, length(ik))]
ik <- ik[-c(1, length(ik))]

frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepSplineB$new(x, df = 11L, intercept = FALSE))$
  plate("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_spline_b(x, deg_free = 11L, complete_set = FALSE, keep_original_cols = TRUE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

expect_equivalent(frec, rec)

frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepSplineB$new(x, df = 11L, intercept = TRUE))$
  plate("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_spline_b(x, deg_free = 11L, complete_set = TRUE, keep_original_cols = TRUE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

expect_equivalent(frec, rec)
