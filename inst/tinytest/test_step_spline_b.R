set.seed(1)

formula <- as.formula(x~y+z)
rows <- 1e5

dat <- data.frame(x = rnorm(rows),
                  y = 1:rows,
                  z = cumsum(rnorm(rows)))
ik <- fquantile(dat$x, probs = seq(0, 1, 0.1))
bk <- ik[c(1, length(ik))]
ik <- ik[-c(1, length(ik))]

frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepSplineB$new(x, df = 11L, intercept = FALSE))$
  prep()$
  bake()$
  data("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_spline_b(x, deg_free = 11L, complete_set = FALSE, keep_original_cols = TRUE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

tinytest::expect_equivalent(frec, rec)

frec = Recipe$new(formula = formula, data = dat)$
  add_step(StepSplineB$new(x, df = 11L, intercept = TRUE))$
  prep()$
  bake()$
  data("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_spline_b(x, deg_free = 11L, complete_set = TRUE, keep_original_cols = TRUE) |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

tinytest::expect_equivalent(frec, rec)
