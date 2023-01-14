test_that("detrend_matrix works", {

  # spec.pgram detrend function
  dt <- function(x) {
    N <- nrow(x)
    t <- 1L:N - (N + 1)/2
    sumt2 <- N * (N^2 - 1)/12
    for (i in 1L:ncol(x))
      x[, i] <- x[, i] - mean(x[, i]) - sum(x[, i] * t) * t/sumt2

    return(x)
  }

  m <- matrix(rnorm(1000), ncol = 10) + 1:100

  expect_equal(detrend_matrix(m), dt(m))

})

test_that("demean_matrix works", {

  m <- matrix(rnorm(1000), ncol = 10) + 1:100

  expect_equal(demean_matrix(m), apply(m, 2, function(x) x - mean(x)))

})

test_that("detrend_and_demean_matrix works", {

  m <- matrix(rnorm(1000), ncol = 10) + 1:100

  expect_equal(detrend_and_demean_matrix(m, TRUE, FALSE),
               detrend_matrix(m))

  expect_equal(detrend_and_demean_matrix(m, FALSE, TRUE),
               demean_matrix(m))

  expect_equal(detrend_and_demean_matrix(m, TRUE, TRUE),
               apply(detrend_matrix(m), 2, function(x) x - mean(x)))

})
