
# no earth tides all methods should be equivalent
data("kennel_2020")

kennel_2020[["et"]] <- 0.0
kennel_2020[["wl"]] <- 0.3 * kennel_2020[["baro"]]


formula <- as.formula(wl~.)
frec = hydrorecipes:::Recipe$new(formula = formula,
                  data = unclass(kennel_2020))$
  add_step(hydrorecipes:::StepBaroHarmonic$new(datetime,
                                wl,
                                baro,
                                et,
                                inverse = FALSE))$
  prep()$
  bake()

expect_equivalent(unlist(frec$get_step_data("barometric_efficiency")),
                            rep(0.3, 4), info = "StepBaroHarmonic works with no Earth tides")


kennel_2020[["wl"]] <- kennel_2020[["wl"]] - kennel_2020[["baro"]]


formula <- as.formula(wl~.)

frec = hydrorecipes:::Recipe$new(formula = formula,
                  data = unclass(kennel_2020))$
  add_step(hydrorecipes:::StepBaroHarmonic$new(datetime,
                                wl,
                                baro,
                                et,
                                inverse = FALSE))$
  prep()$
  bake()
expect_equivalent(unlist(frec$get_step_data("barometric_efficiency")),
                            rep(0.7, 4), info = "StepBaroHarmonic works with no Earth tides")


frec = hydrorecipes:::Recipe$new(formula = formula,
                                 data = unclass(kennel_2020))$
  add_step(hydrorecipes:::StepBaroHarmonic$new(datetime,
                                               wl,
                                               baro,
                                               et,
                                               inverse = TRUE))$
  prep()$
  bake()
expect_equivalent(unlist(frec$get_step_data("barometric_efficiency")),
                            rep(0.3, 4), info = "StepBaroHarmonic works with no Earth tides")



#
# Mod(hydrorecipes:::be_transfer(qM(kennel_2020[, list(wl, baro, et)]), c(3),
#                            TRUE, TRUE, 0.2, 1.0, 1440.0))
#
#
# data("kennel_2020")
#
# # earth tides methods will not be equivalent (ratio should fail)
# kennel_2020[, et := sin(2 * pi / 86400 * 2.0 * as.numeric(datetime))  + 2*sin(2 * pi / 86400 * 1.9324 * as.numeric(datetime)) + rnorm(nrow(kennel_2020), sd = 0.001)]
# # kennel_2020[, et := 0]
# # kennel_2020[, et := et + rnorm(nrow(kennel_2020), sd = 0.0001)]
# kennel_2020[, wl := 0.3 * baro + 0.001 * (sin(2 * pi / 86400 * 2.0 * as.numeric(datetime))  + 2*sin(2 * pi / 86400 * 1.9324 * as.numeric(datetime)))]
#
#
#
# bench::mark({
#
# formula <- as.formula(wl~.)
# frec = Recipe$new(formula = formula,
#                   data = unclass(kennel_2020))$
#   add_step(StepBaroHarmonic$new(datetime,
#                                 wl,
#                                 baro,
#                                 et,
#                                 inverse = FALSE))$
#   prep()$
#   bake()
# unlist(frec$get_step_data("barometric_efficiency"))
# })
#
# Mod(hydrorecipes:::be_transfer(qM(kennel_2020[, list(wl, baro, et)]), 3,
#                            TRUE, TRUE, 0.2, 2.0, 1440.0))
#
#
#
#
# dat <- fread("../../r_scratch/death_valley.csv")[,1:4]
# setnames(dat, c("datetime", "wl", "baro", "et"))
# dat[, datetime := as.POSIXct(datetime, format = "%d/%m/%Y %H:%M", tz = "US/Pacific")]
# Mod(hydrorecipes:::be_transfer(qM(dat[, list(wl, baro, et)]), c(3,3,3),
#                            TRUE, TRUE, 0.4, 0.6, 1440.0 / 5.0))
#
#
# goertzel <- function(signal, k, N) {
#   # signal: Input signal
#   # k: Index of the frequency bin to calculate (0-indexed)
#   # N: Length of the signal
#
#   omega <- 2 * pi * k / N
#   sine <- sin(omega)
#   cosine <- cos(omega)
#   coeff <- 2 * cosine
#
#   Q1 <- 0
#   Q2 <- 0
#
#   for (n in 1:N) {
#     Q0 <- coeff * Q1 - Q2 + signal[n]
#     Q2 <- Q1
#     Q1 <- Q0
#   }
#
#   real_part <- Q1 - Q2 * cosine
#   imag_part <- Q2 * sine
#
#   return(complex(real = real_part, imaginary = imag_part))
# }
#
# plot(butter(4))
# goertzel(kennel_2020$wl, 56.5, nrow(kennel_2020))
# goertzel(kennel_2020$baro, 111, nrow(kennel_2020))
# goertzel(kennel_2020$et, 111, nrow(kennel_2020))
# tmp <- sapply(1:200, function(x) {
#   goertzel(kennel_2020$wl, x, nrow(kennel_2020))/
#     goertzel(kennel_2020$baro, x, nrow(kennel_2020))
#   })
# plot(Mod(tmp), type = 'l')
# tmp1 <- sapply(1:200, function(x) {
#   goertzel(kennel_2020$wl, x, nrow(kennel_2020))
# })
# plot(Mod(tmp1), type = 'l')
#
# tmp2 <- sapply(1:200, function(x) {
#   hydrorecipes:::dft_goertzel(qM(kennel_2020[, list(wl)]), x)
# })
# tmp1 <- sapply(1:150, function(x) {
# Mod(hydrorecipes:::be_dft(qM(kennel_2020[, list(wl, baro, et)]), x))
# })
# plot((tmp1[1,]), type = 'l', log = 'y')
#
# tmp1a <- sapply(1:150, function(x) {
#   hydrorecipes:::be_dft(qM(kennel_2020[1:20000, list(wl, baro, et)]), x)
# })
#
# k <- 1
# test <- sapply(1:200, function(x) {
# hydrorecipes:::be_dft(qM(kennel_2020[(1440 * k):(1440*(8+k)), list(wl, baro, et)]), x)
# })
# plot(Mod(test[1,]), type = 'l')
# plot(Mod(test[2,]), type = 'l', log = 'y')
# abline(v = 8)
# abline(h = 0.001)
#
#
# tmp1b <- sapply(1:150, function(x) {
#   hydrorecipes:::be_dft(qM(kennel_2020[20000:40000, list(wl, baro, et)]), x)
# })
#
# tmp1c <- sapply(1:150, function(x) {
#   hydrorecipes:::be_dft(qM(kennel_2020[40000:60000, list(wl, baro, et)]), x)
# })
#
# tmp1d <- sapply(1:150, function(x) {
#   hydrorecipes:::be_dft(qM(kennel_2020[60000:80000, list(wl, baro, et)]), x)
# })
#
# plot(Mod(tmp1a[2,] + tmp1b[2,] + tmp1c[2,] + tmp1d[2,]), type = 'l')
#
#
# tf <- hydrorecipes:::transfer_pgram(qM(kennel_2020[, list(wl, baro, et)]), spans = c(3), FALSE, FALSE, 0.45)
#
#
# plot(Mod(tf[,1]), type = 'l', col = 'red', xlim = c(0, 200))
# plot(Mod(tf[,2]), type = 'l', col = 'red', xlim = c(0, 200))
#
# plot(Mod(tmp[2,1:200]), type = 'p', log = 'y', xlim = c(0,100))
#
#
#
# tmp <- hydrorecipes:::be_dft(qM(kennel_2020[, list(wl, baro, et)]), 113)
#
# a <- transfer_pgram(qM(kennel_2020[, list(wl, baro, et)]), 3, TRUE, TRUE, 0.1)
#
#
# Mod(dft(qM(kennel_2020[, list(wl)]), 100))/
# Mod(dft(qM(kennel_2020[, list(baro)]), 100))
#
#
# expect_equivalent(unlist(frec$get_step_data("barometric_efficiency"))[2:3],
#                             rep(0.3, 2), tolerance = 0.05, info = "StepBaroHarmonic works small Earth tides")
#
#
# kennel_2020[, et := 2*cos(2 * pi / 1440 * 1.9324 * as.numeric(datetime)) +
#               cos(2 * pi / 1440 * 2.0 * as.numeric(datetime))]
# kennel_2020[, wl := 0.5 * baro + et * 0.001]
# # kennel_2020[, wl := wl - baro]
#
# formula <- as.formula(wl~.)
#
# frec = Recipe$new(formula = formula,
#                   data = unclass(kennel_2020))$
#   add_step(StepBaroHarmonic$new(datetime,
#                                 wl,
#                                 baro,
#                                 et,
#                                 inverse = TRUE))$
#   prep()$
#   bake()
#
# expect_equivalent(unlist(frec$get_step_data("barometric_efficiency"))[2:3],
#                             rep(0.7, 2), tolerance = 0.05, info = "StepBaroHarmonic works small Earth tides")


