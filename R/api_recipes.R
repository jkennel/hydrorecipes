
# new recipe -------------------------------------------------------------------

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' Create a new R6 recipe. This is analogous to the the list structure that the
#' *recipes* package uses.
#'
#' @inheritParams stats::lm
#' @param ... additional arguments to pass to Recipe$new().  This is currently
#' not used.
#'
#' @return R6 recipe object
#' @export
#'
#' @examples
#' rec <- recipe(y~x, data = list(x = rnorm(10), y = rnorm(10)))
#'
recipe <- function(formula, data, ...) {
  Recipe$new(formula, data, ...)
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

# steps ------------------------------------------------------------------------

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_scale
#'
#' @description
#'   Adds a step to scale a data column(s)
#'
#' @param .rec the R6 recipe object.
#' @param terms the unquoted names of the variables to use or a selector
#'   function.  terms replaces the `...` of the recipes package but requires
#'   variables to be included within `c()`.  For example to include variables x
#'   and y you would write `c(x,y)` in the frecipes package.
#' @param role character - the name of the role
#' @param skip logical - should the step be skipped
#' @param na_rm logical - should NA values be removed from calculations
#' @param fun function - the function that is applied to a list or columns of a
#'   data.frame like object.
#' @param n_sd numeric - number of standard deviations for the scaling
#' @param keep_original_cols logical - keep the original columns or replace them
#' @param ... additional arguments
#'
#' @return
#' @export
#'
#' @examples
#'
#' rec <- recipe(y~x, data = list(x = rnorm(10), y = rnorm(10))) |>
#'        step_scale(x)
step_scale <- function(.rec,
                       terms,
                       role = "predictor",
                       skip = FALSE,
                       na_rm = TRUE,
                       fun = collapse::fsd,
                       n_sd = 1L,
                       keep_original_cols = FALSE,
                       ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments_no_rec()
  .rec$add_step(do.call(StepScale$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_center
#' @description
#'   Adds a step to center a data column(s)
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#' rec <- recipe(y~x, data = list(x = rnorm(10), y = rnorm(10))) |>
#'        step_center(x)
step_center <- function(.rec,
                        terms,
                        role = "predictor",
                        skip = FALSE,
                        na_rm = TRUE,
                        fun = collapse::fmean,
                        keep_original_cols = FALSE,
                        ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepCenter$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_normalize
#'
#' @inheritParams step_scale
#'
#' @return
#' @export
#'
#' @examples
#' rec <- recipe(y~x, data = list(x = rnorm(10), y = rnorm(10))) |>
#'        step_normalize(x)
step_normalize <- function(.rec,
                           terms,
                           role = "predictor",
                           skip = FALSE,
                           na_rm = TRUE,
                           keep_original_cols = FALSE,
                           ...){

  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepNormalize$new,
                        env_list))
}


#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_harmonic
#'
#' @description
#'   Add sin and cos terms for harmonic analysis
#'
#' @inheritParams step_scale
#' @param frequency numeric vector - the frequencies of the sin and cos curves
#' @param cycle_size numeric - the period of the sin and cos curves
#' @param starting_value numeric - the starting position of the sin and cos
#'   curves. This may be specified to have more control over the signal phase.
#'
#' @return
#' @export
#'
#' @examples
#' rec <- recipe(y~x, data = list(x = 1:10, y = rnorm(10))) |>
#'        step_harmonic(x,
#'                      frequency = 2.0,
#'                      cycle_size = 4.0,
#'                      starting_value = 0.0)
step_harmonic <- function(.rec,
                          terms,
                          frequency = NA_real_,
                          cycle_size = NA_real_,
                          starting_value = NA_real_,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepHarmonic$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_lead_lag
#'
#' @description
#'   Lag or lead a column or columns.  This requires a sorted and regular time
#'   series.
#'
#' @inheritParams step_scale
#' @param lag integer vector - number of samples to lag or lead. Negative
#'   numbers indicate leading a vector.
#' @param n_shift integer - number of values to shift the starting position when
#'   n_subset is not equal to 0. The value of n_shift has to be less than
#'   `n_subset`.
#' @param n_subset integer - spacing between adjacent samples in the result.
#'
#' @return
#' @export
#'
#' @examples
#' rec <- recipe(y~x, data = list(x = rnorm(10), y = rnorm(10))) |>
#'        step_lead_lag(x, lag = 1)
#' rec <- recipe(y~x, data = list(x = rnorm(10), y = rnorm(10))) |>
#'        step_lead_lag(x, lag = 1, n_subset = 5)
#' rec <- recipe(y~x, data = list(x = rnorm(10), y = rnorm(10))) |>
#'        step_lead_lag(x, lag = 1, n_shift = 2, n_subset = 5)
#'
step_lead_lag <- function(terms,
                          lag,
                          n_shift = 0L,
                          n_subset = 1L,
                          role = "predictor",
                          skip = FALSE,
                          keep_original_cols = FALSE,
                          ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepLeadLag$new,
                        env_list))
}
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#' @title step_dummy
#'
#' @description
#'   dummy encode a factor or factor like variable.
#'
#'
#' @inheritParams step_scale
#' @param one_hot logical - use one hot encoding.
#'
#' @return
#' @export
#'
#' @examples
#' dat <- data.frame(x = qF(sample(1:10, 100, replace = TRUE)),
#'                   y = rnorm(100))
#'
#' rec <- recipe(y~x, data = dat) |>
#'        step_dummy(x, one_hot = FALSE)
#' rec <- recipe(y~x, data = dat) |>
#'        step_dummy(x, one_hot = TRUE)
step_dummy <- function(terms,
                       one_hot = FALSE,
                       role = "predictor",
                       skip = FALSE,
                       keep_original_cols = FALSE,
                       ...) {
  terms <- substitute(terms)
  env_list <- get_function_arguments()
  .rec$add_step(do.call(StepDummy$new,
                        env_list))
}


# prep --------------------------------------------------------------------
prep <- function(.rec, retain = TRUE) {
  .rec$prep(retain)
}


# bake --------------------------------------------------------------------

bake <- function(.rec, new_data = NULL, type = "list") {
  .rec$bake(new_data, type)
}



# formula <- as.formula(y~x)
# data <- data.frame(x = as.numeric(1:10000), y = as.numeric(1:10000))
# dat <- data
# frec4 <- frecipes:::recipe(formula, data) |>
#   step_normalize(x) |>
#   prep() |>
#   bake()
#
# bench::mark(
#   rec1 <- recipes::recipe(formula, data) |>
#     recipes::step_scale(x) |>
#     recipes::prep() |>
#     recipes::bake(new_data = NULL),
#   frec2 = Recipe$new(formula = formula, data = data)$
#     add_step(StepScale$new(x))$
#     prep()$
#     bake(),
#   frec1 <- frecipes:::recipe(formula, data) |>
#     step_scale(x) |>
#     prep() |>
#     bake(),
#   frec3 <- frecipes:::recipe(formula, data) |>
#     step_center(x) |>
#     prep() |>
#     bake(),
#   frec4 <- frecipes:::recipe(formula, data) |>
#     step_normalize(x) |>
#     prep() |>
#     bake(),
#   frec5 <- recipes:::recipe(formula, data) |>
#     recipes::step_normalize(x) |>
#     recipes::prep() |>
#     recipes::bake(new_data = NULL),
#   check = FALSE
#   # relative = TRUE
# )


