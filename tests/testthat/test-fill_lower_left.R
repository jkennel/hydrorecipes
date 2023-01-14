# test_that("fill_lower_left works", {
#
#   m <- matrix(rnorm(900), ncol = 9)
#   m[,4] <- 0.0
#   m[,7] <- 0.0
#   m[,8] <- 0.0
#   a <- fft_matrix(m, n_new = nrow(m))
#
#   b <- fill_lower_left(a, n_col = 3, start = 0)
#
#   expect_equal(a[,2], Conj(b[,4]))
#   expect_equal(a[,3], Conj(b[,7]))
#   expect_equal(a[,6], Conj(b[,8]))
#
# })
