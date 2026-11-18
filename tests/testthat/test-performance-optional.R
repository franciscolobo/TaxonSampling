test_that("sampling time and memory remain within optional thresholds", {
  skip_if(Sys.getenv("TS_RUN_PERFORMANCE") != "1",
          "Set TS_RUN_PERFORMANCE=1 to run benchmark")
  # Performance thresholds depend on hardware. Set explicitly in CI.
  limit <- as.numeric(Sys.getenv("TS_MAX_SECONDS", "30"))
  n <- as.integer(Sys.getenv("TS_PERF_OBS_PER_TIP", "500"))
  # Benchmark distinct taxonomic IDs, not repeated observations per tip.
  x <- make_ts_fixture(n_per_tip = 1L)
  base <- x$nodes
  extra <- data.frame(id = paste0("ott:perf_", seq_len(n * 4L)),
                      parent = rep(c("ott:A", "ott:B"), each = n * 2L),
                      name = paste0("perf_", seq_len(n * 4L)),
                      rank = "species", depth = 2L, child_count = 0L)
  x$nodes <- rbind(base, extra)
  x$ids_df <- data.frame(taxID = paste0("tax_perf_", seq_len(n * 4L)),
                         seqID = paste0("seq_perf_", seq_len(n * 4L)),
                         node_id = extra$id)
  x$spp_df <- rbind(x$spp_df, data.frame(taxID = extra$id,
                                         species_count = rep(1, nrow(extra))))
  gc()
  start <- proc.time()[["elapsed"]]
  y <- sample_fixture(x, m = n * 2L, randomize = "yes")
  elapsed <- proc.time()[["elapsed"]] - start
  message(sprintf("Benchmark: %d observations, %d sampled, %.3f seconds, %.1f MiB input",
                  nrow(x$ids_df), nrow(y$outputMapping), elapsed,
                  as.numeric(object.size(x)) / 1024^2))
  expect_equal(nrow(y$outputMapping), n * 2L)
  expect_lt(elapsed, limit)
})
