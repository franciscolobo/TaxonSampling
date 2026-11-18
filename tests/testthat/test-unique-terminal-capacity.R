test_that("terminal taxa have capacity one without replacement", {
  x <- make_ts_fixture(n_per_tip = 7L)
  for (method in c("diversity", "balanced")) {
    for (randomize in c("no", "yes", "after_first_round")) {
      if (method == "balanced" && randomize == "after_first_round") next
      for (sampling in c("agnostic", "known_species")) {
        set.seed(2026)
        y <- sample_fixture(x, m = 4L, method = method,
                            randomize = randomize, sampling = sampling)
        expect_equal(length(y$outputNodes), 4L)
        expect_equal(anyDuplicated(y$outputNodes), 0L)
        expect_equal(anyDuplicated(y$outputIDs), 0L)
        expect_error(sample_fixture(x, m = 5L, method = method,
                                    randomize = randomize, sampling = sampling),
                     "unique eligible taxonomic nodes")
      }
    }
  }
})

test_that("duplicate input records never increase terminal capacity", {
  x <- make_ts_fixture(n_per_tip = 2L)
  x$ids_df <- rbind(x$ids_df, x$ids_df, x$ids_df)
  y <- sample_fixture(x, m = 4L)
  expect_equal(length(unique(y$outputNodes)), 4L)
  expect_error(sample_fixture(x, m = 5L), "unique eligible taxonomic nodes")
})

test_that("terminal taxa can repeat only with explicit replacement", {
  x <- make_ts_fixture(n_per_tip = 5L)
  y <- sample_fixture(x, m = 12L, replacement = TRUE, randomize = "yes")
  expect_equal(length(y$outputNodes), 12L)
  expect_lt(length(unique(y$outputNodes)), 12L)
})
