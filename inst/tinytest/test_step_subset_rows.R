dat <- data.frame(x = as.numeric(1:200),
                  y = rep(0.01, 200))
formula <- as.formula(y~x)


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
# frec1 = recipe(formula = formula, data = dat) |>
#   step_subset_rows(x, row_numbers = c(1, 5, 10)) |>
#   plate("dt")


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
