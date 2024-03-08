
# test calculation -------------------------------------------------------------
ac_1 <-  be_acworth_calc_cpp(s2_at = 7.461,
                             s2_et = 224.640,
                             s2_gw = 4.086,
                             m2_gw = 0.471,
                             m2_et = 492.526,
                             d_phase = -56.709 * pi / 180,
                             inverse = TRUE)
ac_2 <- be_acworth_calc_cpp(s2_at = 6.164,
                            s2_et = 270.463,
                            s2_gw = 0.329,
                            m2_gw = 0.225,
                            m2_et = 551.572,
                            d_phase = -71.726 * pi / 180,
                            inverse = TRUE)
ac_3 <- be_acworth_calc_cpp(s2_at = 5.897,
                            s2_et = 234.478,
                            s2_gw = 5.536,
                            m2_gw = 0.773,
                            m2_et = 558.075,
                            d_phase = -70.393 * pi / 180,
                            inverse = TRUE)
tinytest::expect_equivalent(ac_1, 0.563,  tolerance = 1e-3,
                            info = "be_acworth_cpp works")
tinytest::expect_equivalent(ac_2, 0.059,  tolerance = 1e-3,
                            info = "be_acworth_cpp works")


# test peaks -------------------------------------------------------------------

p <- frecipes:::get_peaks(as.numeric(1:1000), 50.87, 70.1)

tinytest::expect_equivalent(p, c(51 - 1, 70 - 1),
                            info = "get_peaks works")


# test inverse -----------------------------------------------------------------
n <- 1000
x <- rnorm(n)
y <- -(x - x * 0.4)
z <- rep(0.0, n)

dat <- matrix(c(y, x, z), ncol = 3)
be_in <- be_acworth_cpp(dat, 3, FALSE, FALSE, 0.1,
                     TRUE, 1.93/86400, 2/86400, 86400)

y <- x * 0.4
dat <- matrix(c(y, x, z), ncol = 3)
be_le <- be_acworth_cpp(dat, 3, FALSE, FALSE, 0.1,
                     FALSE, 1.93/86400, 2/86400, 86400)

tinytest::expect_equivalent(p, c(51 - 1, 70 - 1),
                            info = "inverse parameter works")
