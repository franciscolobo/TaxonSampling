test_that("duplicates and multiple IDs per terminal node never inflate capacity", {
  x <- make_ts_fixture(n_per_tip = 5L)
  for (method in c("diversity", "balanced")) {
    for (randomize in c("no", "yes", "after_first_round")) {
      if (method == "balanced" && randomize == "after_first_round") next
      for (sampling in c("agnostic", "known_species")) {
        set.seed(14)
        y <- sample_fixture(x, m = 4L, method = method,
                            randomize = randomize, sampling = sampling)
        expect_length(unique(y$outputNodes), 4L)
        expect_length(unique(y$outputIDs), 4L)
        expect_equal(nrow(y$outputMapping), 4L)
      }
    }
  }
  expect_error(sample_fixture(x, m = 5L), "unique eligible taxonomic nodes")
  x$ids_df <- rbind(x$ids_df, x$ids_df)
  expect_length(unique(sample_fixture(x, m = 4L)$outputNodes), 4L)
})

test_that("direct internal observations count once and exclusions reduce capacity", {
  x <- make_ts_fixture(n_per_tip = 4L, direct_root = TRUE)
  expect_length(unique(sample_fixture(x, m = 5L)$outputNodes), 5L)
  expect_error(sample_fixture(x, m = 6L), "unique eligible")
  expect_error(sample_fixture(x, m = 4L, ignoreIDs = "ott:A"),
               "unique eligible")
})

test_that("replacement explicitly permits repeated native nodes", {
  x <- make_ts_fixture(n_per_tip = 3L)
  y <- sample_fixture(x, m = 15L, replacement = TRUE, randomize = "yes")
  expect_length(y$outputNodes, 15L)
  expect_true(length(unique(y$outputNodes)) <= 4L)
})

test_that("ambiguous input TaxID mapping is rejected", {
  x <- make_ts_fixture(n_per_tip = 1L)
  x$ids_df <- rbind(x$ids_df, transform(x$ids_df[1L, ], node_id = "ott:B1"))
  expect_error(sample_fixture(x, m = 2L), "multiple native nodes")
})
