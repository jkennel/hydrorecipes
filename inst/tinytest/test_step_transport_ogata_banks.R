
R <- 1.0
D <- 0.1
k <- 0.0 / R
v <- 1.0
c0 <- 1.0
x <- 1.0
t <- 1.4

ob_1 <- frecipes:::ogata_banks_ind(D, v, c0, x, t)

ob_2 <- frecipes:::ogata_banks_decay_ind(c0 = c0,
                                       v = v,
                                       D = D,
                                       R = R,
                                       k = k,
                                       x = x,
                                       t = t)


tinytest::expect_equivalent(ob_1, ob_2,
                            info = "Ogata Banks without decay and retardation are equal")




# value from https://www.civil.uwaterloo.ca/jrcraig/pdf/OgataBanks.xlsm
tinytest::expect_equivalent(0.838457815, ob_2,
                            tolerance = 1e-4,
                            info = "Ogata Banks without decay and retardation are equal")


R <- 1.25
D <- 0.1
k <- 0.2 / R
v <- 1.0
c0 <- 1.0
x <- 1.0
t <- 1.4

ob_2 <- frecipes:::ogata_banks_decay_ind(c0 = c0,
                                         v = v,
                                         D = D,
                                         R = R,
                                         k = k,
                                         x = x,
                                         t = t)

# value from https://www.civil.uwaterloo.ca/jrcraig/pdf/OgataBanks.xlsm
tinytest::expect_equivalent(0.587007929, ob_2,
                            tolerance = 1e-4,
                            info = "Ogata Banks with decay and retardation are equal")

R <- 1
D <- 0.1
k <- 0.0 / R
v <- 1.0
c0 <- 1.0
x <- c(-1.0, 0.0, 10000.0, 0.0, 10000.0)
t <- c(-1.0, 0.0, 0.0, 10000.0, 10000.0)

ob_2 <- frecipes:::ogata_banks_decay_vec(c0 = c0,
                                         v = v,
                                         D = D,
                                         R = R,
                                         k = k,
                                         x = x,
                                         t = t)

tinytest::expect_equivalent(c(0.0, 0.0, 0.0, 1.0, 0.0), ob_2,
                            info = "Ogata Banks with short and long values")






# R <- 1.25
# D <- 0.1
# k <- 0.2 / R
# v <- 1.0
# c0 <- 1.0
# x <- 10000.0
# t <- 1.4
#
# ob_2 <- frecipes:::ogata_banks_decay_ind(c0 = c0,
#                                          v = v,
#                                          D = D,
#                                          R = R,
#                                          k = k,
#                                          x = x,
#                                          t = t)
#
# bench::mark(
#
# ob_2 <- ogata_banks_decay_vec(c0 = c0,
#                       v = v,
#                       D = D,
#                       R = R,
#                       k = k,
#                       x = as.numeric(1:1000000)/1000000.0,
#                       t = rep(1.0, 1000000))
# )
# ob_2
