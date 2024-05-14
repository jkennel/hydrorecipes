set.seed(123)

formula <- as.formula(x~y+z)
rows <- 11

dat <- data.frame(x = rep(1, rows),
                  y = 1:rows,
                  z = cumsum(rnorm(rows)))

frec = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepKernelFilter$new(z,
                                           kernel = list(rep(1, 3)/3),
                                           align = "center"))$
  plate("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
  recipes::step_window(z, size = 3, statistic = "mean") |>
  recipes::prep() |>
  recipes::bake(new_data = NULL)

expect_equivalent(frec[2:10, 4]$kernel_filter_z,
                  rec[2:10, "z"]$z)

