test_that("index_from_i_j works", {
  expect_equal(index_from_i_j(2, 3, 5), 13)


  expect_equal(index_from_i_j(2, 3, 4),
               index_from_j_i(3, 2, 4))

  #  0,  1,  2,    3,  4,
  #  5,  6,  7,    8,  9,
  # 10, 11, 12, *13*, 14,
  # 15, 16, 17,   18, 19,
  # 20, 21, 22,   23, 24

  expect_error(index_from_i_j(5, 3, 4))
  expect_error(index_from_i_j(5, 5, 4))
  expect_error(index_from_i_j(2, 5, 4))

})

test_that("index_from_i_j works", {

  expect_equal(index_from_j_i(2, 3, 5), 17)
  #  0,  1,    2,  3,  4,
  #  5,  6,    7,  8,  9,
  # 10, 11,   12, 13, 14,
  # 15, 16, *17*, 18, 19,
  # 20, 21,   22, 23, 24

  expect_error(index_from_j_i(5, 3, 4))
  expect_error(index_from_j_i(5, 5, 4))
  expect_error(index_from_j_i(2, 5, 4))
})
