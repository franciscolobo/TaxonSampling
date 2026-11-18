test_that("allocator supports positional unnamed vectors", {
  allocate <- getFromNamespace(".ts_allocate", "TaxonSampling")
  expect_error(
    allocate(6L, c(0L, 0L), c(2L, 2L), c(1, 1), "diversity", "no"),
    "capacity"
  )
  ans <- allocate(2L, c(0L, 0L), c(2L, 2L), c(1, 1), "diversity", "no")
  expect_equal(unname(ans), c(1L, 1L))
  expect_equal(sum(ans), 2L)
})

test_that("allocator rejects mixed named and unnamed vectors", {
  allocate <- getFromNamespace(".ts_allocate", "TaxonSampling")
  expect_error(
    allocate(1L, c(A = 0L, B = 0L), c(2L, 2L), c(1, 1),
             "diversity", "no"),
    "all be named or all be unnamed"
  )
})
