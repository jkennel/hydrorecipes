test_that("modified_daniell works", {

  n <- 100L
  a <- modified_daniell(c(3,3,7), n)
  b <- kernel(coef = "modified.daniell", c(3,3,7))
  b <- c(b[0L:b$m], rep_len(0,n-2L*b$m-1L), b[-b$m:-1L])
  expect_equal(a,b)


  a <- modified_daniell(3, n)
  b <- kernel(coef = "modified.daniell", 3)
  b <- c(b[0L:b$m], rep_len(0,n-2L*b$m-1L), b[-b$m:-1L])
  expect_equal(a,b)

})
