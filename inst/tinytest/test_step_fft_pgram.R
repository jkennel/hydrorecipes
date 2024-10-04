formula <- as.formula(y~.)

dat <- data.frame(x = rnorm(200),
                  y = rnorm(200))

formula <- as.formula(.~x+y)
frec1 = hydrorecipes:::Recipe$new(formula = formula, data = dat)$
  add_step(hydrorecipes:::StepPgram$new(c(x, y)))$
  plate("df")

frec2 = recipe(formula = formula, data = dat) |>
  step_fft_pgram(c(x,y)) |>
  plate("df")

expect_equivalent(frec1, frec2,
                            info = "R6 and hydrorecipes api are equivalent")


n <- 10000
x <- rnorm(n)
y <- x * 0.2
m <- matrix(c(y,x), ncol = 2)
gain <- as.numeric(Mod(hydrorecipes:::transfer_pgram(m,
                                                 spans = 3,
                                                 detrend = FALSE,
                                                 demean = FALSE,
                                                 taper = 0.1)))

expect_equivalent(gain, rep(0.2, n),
                            info = "transfer_pgram gives the right gain")


n_groups <- 50
gain = as.numeric(Mod(hydrorecipes:::transfer_pgram_smooth(m,
                                                       spans = 3,
                                                       detrend = FALSE,
                                                       demean = FALSE,
                                                       taper = 0.1,
                                                       # power = 3,
                                                       n_groups = n_groups)))
expect_equivalent(gain, rep(0.2, n_groups),
                            info = "transfer_pgram_smooth gives the right gain")



data(kennel_2020)
formula <- as.formula(datetime~baro+wl+et)
formula_1 <- as.formula(wl~baro + et)

frec2 = recipe(formula = formula, data = kennel_2020) |>
  step_fft_transfer_welch(formula = formula_1, length_subset = 1000,
                          window = hydrorecipes:::window_rectangle(1000),
                          time_step = 60) |>
  step_fft_transfer_pgram(formula = formula_1, spans = c(3, 3), time_step = 60) |>
  step_fft_transfer_experimental(formula = formula_1, n_groups = 600,
                                 spans = c(3,3), taper = 0.05, time_step = 60) |>
  prep("df") |>
  bake()

frec2$get_elapsed_times()
tf <- frec2$get_transfer_data(type = "dt")
expect_equivalent(class(tf),
                  c("data.table", "data.frame"),
                  info = "extracting the transfer function returns a data.table")

expect_equivalent(nrow(tf), 167040L,
                  info = "transfer number of rows is correct")


tf <- frec2$get_transfer_data(type = "df")
expect_equivalent(class(tf),
                  c("data.frame"),
                  info = "extracting the transfer function returns a data.frame")

#
#
#
# library(ggplot2)
# ggplot(a, aes(x = frequency, y = Mod(value), color = id, group = id)) +
#   geom_line(alpha = 0.5) +
#   # geom_point() +
#   facet_wrap(variable~., ncol = 1) +
#   scale_x_log10() +
#   scale_y_continuous(limits = c(0, 1))


# plot(Mod(fft_transfer_experimental_1)~frequency, a, type = "l", log = "x", ylim = c(0, 1))
# abline(v = 2)
# abline(v = 1)

#
# frec2 = recipe(formula = formula, data = kennel_2020) |>
#   step_fft_pgram(c(wl, baro, et), spans = c(3), time_step = 60) |>
#   prep("df") |>
#   bake()
# a <- collapse::qDT(frec2$get_step_data(field_name = "fft_result")[[1]])
#

