set.seed(1)

formula <- as.formula(x~y+z)
rows <- 1e5

dat <- data.frame(x = rep(1, rows),
                  y = 1:rows,
                  z = cumsum(rnorm(rows)))

frec = Recipe$new(formula = formula, data = dat)$
    add_step(StepKernelFilter$new(z,
                                  kernel = list(rep(1, 1001)/1001), align = "center"))$
    plate("tbl")

rec  = recipes::recipe(formula = formula, data = dat) |>
    recipes::step_window(z, size = 1001, statistic = "mean",) |>
    recipes::prep() |>
    recipes::bake(new_data = NULL)

tinytest::expect_equivalent(frec[10000:20000,4], rec[10000:20000,'z'])

