test_that("allocation conserves total and respects lower/capacity", {
  # Internal allocator is tested separately from the recursive engine.
  allocate <- getFromNamespace(".ts_allocate", "TaxonSampling")
  for (method in c("diversity", "balanced")) {
    for (mode in c("no", "yes", "after_first_round")) {
      if (method == "balanced" && mode == "after_first_round") next
      for (m in 3:12) {
        set.seed(m)
        a <- allocate(m = m, lower = c(a = 1L, b = 0L, c = 2L),
                      capacity = c(a = 4L, b = 5L, c = 5L),
                      weights = c(a = 1, b = 3, c = 2),
                      method = method, randomize = mode)
        expect_equal(sum(a), m)
        expect_true(all(a >= c(a = 1, b = 0, c = 2)))
        expect_true(all(a <= c(a = 4, b = 5, c = 5)))
      }
    }
  }
})

test_that("known-species weighting changes stochastic selection probability", {
  # One draw, two clades. Probabilities should approximate 0.8 and 0.2.
  allocate <- getFromNamespace(".ts_allocate", "TaxonSampling")
  set.seed(999)
  n <- 2500L
  chosen <- replicate(n, {
    a <- allocate(1L, c(A = 0L, B = 0L), c(A = 1L, B = 1L),
                  c(A = 8, B = 2), "diversity", "yes")
    a[["A"]] == 1L
  })
  p <- mean(chosen)
  # Broad tolerance prevents flaky tests while detecting uniform allocation.
  expect_gt(p, 0.76)
  expect_lt(p, 0.84)
})

test_that("agnostic stochastic selection is approximately symmetric", {
  allocate <- getFromNamespace(".ts_allocate", "TaxonSampling")
  set.seed(1001)
  chosen <- replicate(2500L, {
    a <- allocate(1L, c(A = 0L, B = 0L), c(A = 1L, B = 1L),
                  c(A = 1, B = 1), "diversity", "yes")
    a[["A"]] == 1L
  })
  expect_lt(abs(mean(chosen) - 0.5), 0.04)
})

test_that("impossible allocator requests are rejected", {
  allocate <- getFromNamespace(".ts_allocate", "TaxonSampling")
  expect_error(allocate(1L, c(1L, 1L), c(2L, 2L), c(1, 1),
                        "diversity", "no"), "Required")
  expect_error(allocate(6L, c(0L, 0L), c(2L, 2L), c(1, 1),
                        "diversity", "no"), "capacity")
})
