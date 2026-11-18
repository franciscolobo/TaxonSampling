test_that("sampling with replacement permits m above unique observations", {
  # Verify that repeated observations are allowed when replacement = TRUE.
  x <- make_ts_fixture(n_per_tip = 1L)
  y <- sample_fixture(x, m = 12L, replacement = TRUE, randomize = "yes")
  expect_equal(nrow(y$outputMapping), 12L)
  expect_true(all(y$outputMapping$taxID %in% x$ids_df$taxID))
})

test_that("required ancestor with direct observation is not forced", {
  # Statistical regression: A is required, but all 3 available descendants
  # should be eligible. Current implementation forces the direct A node.
  x <- make_ts_fixture(n_per_tip = 1L)
  x$ids_df <- rbind(x$ids_df,
    data.frame(taxID = "direct_A", seqID = "direct_A_seq", node_id = "ott:A"))
  set.seed(1)
  draws <- replicate(80L, {
    y <- sample_fixture(x, m = 1L, requireIDs = "ott:A", randomize = "yes")
    y$outputMapping$node_id[[1L]]
  })
  expect_true(any(draws != "ott:A"))
})
