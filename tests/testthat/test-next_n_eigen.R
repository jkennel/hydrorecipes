test_that("next_n_eigen works", {

  n <- 100
  expect_equal(stats::nextn(n), next_n_eigen(n))
  n <- 79
  expect_equal(stats::nextn(n), next_n_eigen(n))
  n <- 2^10
  expect_equal(stats::nextn(n), next_n_eigen(n))

})
