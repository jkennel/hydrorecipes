
test_that("spec_taper works", {

  a <- spec_taper(100, 0.1)
  b <- spec.taper(rep(1, 100), 0.1)
  expect_equal(a, b)

  expect_error(spec_taper(100, 0.6))
  expect_error(spec_taper(100, -0.01))

  # bench::mark(spec_taper(1e5, 0.1),
  #             spec.taper(rep(1, 1e5), 0.1))

})
