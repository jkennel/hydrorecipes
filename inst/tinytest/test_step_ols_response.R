data("kennel_2020")
kennel_2020[, datetime := as.numeric(datetime)]

formula <- as.formula(wl~.)
n_knots <- 12
deg_free <- 27
max_lag <- 1 + 720

frec = Recipe$new(formula = formula, data = unclass(kennel_2020))$
  add_step(StepDistributedLag$new(baro,
                                  knots = frecipes:::log_lags_arma(n_knots, max_lag)))$
  add_step(StepSplineB$new(datetime, df = deg_free, intercept = FALSE))$
  add_step(StepIntercept$new())$
  add_step(StepDropColumns$new(baro))$
  add_step(StepDropColumns$new(datetime))$
  add_step(StepOlsResponse$new(formula))$
  prep()$
  bake()

resp <- frec$get_response_data('dt')[variable == "cumulative"]



tinytest::expect_equivalent(resp[.N]$value, 0.879, tolerance = 1e-2)

data("kennel_2020")
kennel_2020[, datetime := as.numeric(datetime)]

formula <- as.formula(wl~.)
n_knots <- 100
deg_free <- 27
max_lag <- 1 + 720

frec = Recipe$new(formula = formula, data = unclass(kennel_2020))$
  add_step(StepLeadLag$new(baro, lag = frecipes:::log_lags_arma(n_knots, max_lag), n_shift = 0, n_subset = 1))$
  add_step(StepSplineB$new(datetime, df = deg_free, intercept = FALSE))$
  add_step(StepIntercept$new())$
  add_step(StepDropColumns$new(baro))$
  add_step(StepDropColumns$new(datetime))$
  add_step(StepOlsResponse$new(formula))$
  prep()$
  bake()

resp <- frec$get_response_data('dt')[variable == "cumulative"]
resp




frec = Recipe$new(formula = formula, data = unclass(kennel_2020))$
  add_step(StepDistributedLag$new(c(baro),
                                  knots = frecipes:::log_lags_arma(n_knots, max_lag)))$
  # add_step(StepSplineB$new(datetime, df = deg_free, intercept = FALSE))$
  # add_step(StepIntercept$new())$
  # add_step(StepDropColumns$new(baro))$
  # add_step(StepDropColumns$new(datetime))$
  add_step(StepOlsResponse$new("wl"))$
  prep()$
  bake()

resp <- frec$get_response_data('dt')[variable == "cumulative"]


plot(value~x, resp, type = 'l', col = 'red', log = 'x')


# bench::mark(
#
#   {rec <- recipes::recipe(formula, kennel_2020) |>
#     hydrorecipes::step_distributed_lag(baro, knots = hydrorecipes::log_lags(n_knots, max_lag)) |>
#     recipes::step_spline_b(datetime, deg_free = deg_free) |>
#     # recipes::step_rm(baro) |>
#     recipes::prep() |>
#     recipes::bake(new_data = NULL);
#
#   fit <- lm(wl~., rec)}
# )

# library(rsk)
# library(data.table)
# library(recipes)
# library(hydrorecipes)
#
# rsk_wl <- function(files) {
#
#   start <- as.POSIXct('2016-08-18', tz = 'UTC')
#   end   <- as.POSIXct('2016-10-13 12:00:00', tz = 'UTC')
#
#   dat <- list()
#   for (i in seq_along(files)) {
#
#     fn <- files[i]
#     d <- Rsk$new(fn)$data
#     d <- rsk::rename_data(d)
#     d <- d[as.numeric(datetime) %% 1 == 0]
#     d <- d[, list(datetime, pressure)]
#     d <- d[datetime %between% c(start, end)]
#
#     col_name <- gsub("_wl.rsk", "", basename(fn))
#     col_name <- gsub(".rsk", "", basename(fn))
#
#     d[, file_name := col_name]
#
#     dat[[i]] <- d
#
#   }
#
#   dat <- rbindlist(dat)
#   dcast(dat, datetime~file_name, value.var = "pressure")
#
# }
#
#
# path <- "../../r_scratch/"
# fn <- file.path(path, c("rd_130_wl.rsk", "rd_130_baro.rsk"))
#
# kennel_2020 <- rsk_wl(fn)
# setnames(kennel_2020, c("datetime", "baro", "wl"))
# start <- as.POSIXct('2016-08-18', tz = 'UTC')
# end   <- as.POSIXct('2016-10-13 12:00:00', tz = 'UTC')
# kennel_2020 <- kennel_2020[datetime %between% c(start, end)]
# formula <- as.formula(wl~.)
# n_knots <- 15
# deg_free <- 27
# max_lag <- 1 + 86400
#
# kennel_2020[, datetime := as.numeric(datetime)]
# rec <- recipes::recipe(formula, kennel_2020) |>
#   hydrorecipes::step_distributed_lag(baro, knots = hydrorecipes::log_lags(n_knots, max_lag)) |>
#   recipes::step_spline_b(datetime, deg_free = deg_free) |>
#   # recipes::step_rm(baro) |>
#   recipes::prep() |>
#   recipes::bake(new_data = NULL)
#
#
# fit <- lm(wl~., rec)
#
# summary(fit)
# (fit$coefficients)
#
# frec = Recipe$new(formula = formula, data = unclass(kennel_2020))$
#   add_step(StepDistributedLag$new(baro,
#                                   knots = frecipes:::log_lags_arma(n_knots, max_lag)))$
#   add_step(StepSplineB$new(datetime, df = deg_free, intercept = FALSE))$
#   add_step(StepIntercept$new())$
#   add_step(StepDropColumns$new(baro))$
#   add_step(StepDropColumns$new(datetime))$
#   add_step(StepOlsResponse$new(wl))$
#   prep()$
#   bake()
#
#
#
# aa <- qDT(frec$steps[[7]]$response_data)
# a <- aa[variable == "cumulative"]
# plot(value~x, a, type = 'l', col = 'red')
#
# #
# # tinytest::expect_equivalent(frec, rec)
# #
# # plot(wl~datetime, kennel_2020, type = 'l')
# # plot(baro~datetime, kennel_2020, type = 'l')
# #
#
# tmp <- qDT(frec$result)
# tmp[, datetime := NULL]
# tmp[, baro := NULL]
# tmp[, intercept := NULL]
# fit2 <- lm(wl~., tmp)
# summary(fit2)
# (fit2$coefficients)
#
#
# knots <- as.numeric(log_lags_arma(n_knots, max_lag))
# rng = 0:(max_lag)
# one_n = c(1L, n_knots)
#
# basis_matrix <- qM(n_spline_list(rng, 0L, 3L, knots[-one_n],
#                                    knots[one_n], TRUE, FALSE,
#                                    0L, FALSE))
# co <- fit2$coefficients
# n <- length(co)
# co <- co[-c(1, (n-deg_free + 1):n)]
# plot(cumsum(basis_matrix %*% co), type ='l')
#
#
#
#
# tmp <- data.table(x = 0:1e5, y = cumsum(rnorm(1e5)))
# tmp2 <- cbind(x = 0:1e5, )
# summary(f1 <- lm(y~splines::ns(x, knots = c(1,100, 1000, 50000), Boundary.knots = c(0, 1e5), intercept = TRUE)-1, tmp))
# summary(f2 <- lm(y~splines2::naturalSpline(x, knots = c(1,100, 1000, 50000), Boundary.knots = c(0, 1e5), intercept = TRUE)-1, tmp))
#
# splines::ns(x, knots = c(1,100, 1000, 50000), Boundary.knots = c(0, 1e5), intercept = TRUE) %*% coefficients(f1)
# splines2::naturalSpline(x, knots = c(1,100, 1000, 50000), Boundary.knots = c(0, 1e5), intercept = TRUE) %*% coefficients(f2)

