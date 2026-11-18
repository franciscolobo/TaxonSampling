test_that("exact size and mapping integrity across available strategies", {
  x <- make_ts_fixture()
  for (method in c("diversity", "balanced")) {
    for (randomize in c("no", "yes", "after_first_round")) {
      if (method == "balanced" && randomize == "after_first_round") next
      for (sampling in c("agnostic", "known_species")) {
        set.seed(42)
        y <- sample_fixture(x, m = 4L, method = method,
                            randomize = randomize, sampling = sampling)
        expect_equal(nrow(y$outputMapping), 4L)
        expect_identical(y$outputMapping$sample_index, seq_len(4L))
        expect_equal(length(unique(y$outputMapping$taxID)), 4L)
        expect_equal(length(unique(y$outputMapping$seqID)), 4L)
        expect_true(all(y$outputMapping$taxID %in% x$ids_df$taxID))
        expect_true(all(y$outputMapping$node_id %in% x$nodes$id))
        expect_true(all(paste(y$outputMapping$taxID, y$outputMapping$node_id) %in%
                        paste(x$ids_df$taxID, x$ids_df$node_id)))
      }
    }
  }
})

test_that("capacity and oversubscription are handled", {
  x <- make_ts_fixture(n_per_tip = 2L)
  y <- sample_fixture(x, m = 4L)
  expect_equal(nrow(y$outputMapping), 4L)
  expect_equal(length(unique(y$outputMapping$taxID)), 4L)
  expect_error(sample_fixture(x, m = 20L), "unique eligible")
  expect_equal(as.integer(table(y$outputMapping$node_id)), rep(1L, 4L))
})

test_that("diversity first-round covers all available root children", {
  x <- make_ts_fixture(n_per_tip = 4L)
  for (mode in c("no", "after_first_round")) {
    for (seed in 1:10) {
      set.seed(seed)
      y <- sample_fixture(x, m = 2L, method = "diversity", randomize = mode)
      expect_true(covered(y, "ott:A"))
      expect_true(covered(y, "ott:B"))
    }
  }
})

test_that("ignoreIDs excludes entire clades", {
  x <- make_ts_fixture()
  for (method in c("diversity", "balanced")) {
    y <- sample_fixture(x, m = 2L, method = method,
                        ignoreIDs = "ott:A", randomize = "yes")
    expect_false(covered(y, "ott:A"))
    expect_true(all(y$outputMapping$node_id %in% c("ott:B1", "ott:B2")))
  }
  expect_error(sample_fixture(x, m = 1L, ignoreIDs = "ott:unknown"), "absent")
})

test_that("required disjoint clades are represented", {
  x <- make_ts_fixture()
  for (seed in 1:10) {
    set.seed(seed)
    y <- sample_fixture(x, m = 2L, randomize = "yes",
                        requireIDs = c("ott:A1", "ott:B2"))
    expect_true(covered(y, "ott:A1"))
    expect_true(covered(y, "ott:B2"))
  }
  expect_error(sample_fixture(x, m = 1L,
    requireIDs = c("ott:A1", "ott:B2")), "smaller")
  expect_error(sample_fixture(x, m = 2L, requireIDs = "ott:A1",
    ignoreIDs = "ott:A"), "excluded")
})

test_that("ancestor requirement can be fulfilled by a descendant", {
  # Regression test: currently expected to reveal direct-node requirement bug.
  x <- make_ts_fixture(direct_root = TRUE)
  y <- sample_fixture(x, m = 1L, requireIDs = "ott:root",
                      ignoreIDs = NULL, randomize = "no")
  expect_true(covered(y, "ott:root"))
  # This test uses a direct root observation and a root requirement; it cannot
  # demand a particular selected taxon. The stronger regression is below.
})

test_that("requirements on internal nodes do not force direct observation", {
  # A required internal clade can be represented by any descendant.
  # A direct observation at A must not be mandatory when B is also required.
  x <- make_ts_fixture(direct_root = FALSE)
  x$ids_df <- rbind(x$ids_df, data.frame(
    taxID = "direct_A", seqID = "direct_A_seq", node_id = "ott:A"))
  # m=2: one from A and one from B; the A observation need not be direct.
  # With an explicitly required descendant, direct_A must not be selected.
  y <- sample_fixture(x, m = 2L, randomize = "no",
                      requireIDs = c("ott:A", "ott:A1", "ott:B1"))
  expect_true(covered(y, "ott:A1"))
  expect_true(covered(y, "ott:B1"))
  expect_false("direct_A" %in% y$outputIDs)
})

test_that("seed reproducibility holds for mapping and metadata", {
  x <- make_ts_fixture()
  set.seed(2026)
  a <- sample_fixture(x, m = 4L, randomize = "yes")
  set.seed(2026)
  b <- sample_fixture(x, m = 4L, randomize = "yes")
  expect_identical(a$outputMapping, b$outputMapping)
  expect_identical(a$ts.params, b$ts.params)
})

test_that("arbitrary backbone IDs and input taxIDs are not conflated", {
  x <- make_ts_fixture(prefix = "synthetic:node_")
  y <- sample_fixture(x, m = 4L)
  expect_true(all(startsWith(y$outputMapping$node_id, "synthetic:node_")))
  expect_true(all(startsWith(y$outputMapping$taxID, "tax_")))
  expect_false(any(y$outputMapping$taxID %in% x$nodes$id))
})

test_that("invalid and infeasible arguments fail clearly", {
  x <- make_ts_fixture()
  expect_error(sample_fixture(x, m = 0L), "positive integer")
  expect_error(sample_fixture(x, m = 1.5), "positive integer")
  expect_error(sample_fixture(x, m = 2L, method = "balanced",
                              randomize = "after_first_round"), "only for diversity")
  expect_error(sample_fixture(x, m = 2L, requireIDs = "ott:absent"), "absent")
  expect_error(TaxonSampling::run_TS(x, "unknown", 2L, verbose = FALSE), "not an observed")
})
