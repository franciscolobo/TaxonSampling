test_that("all ten repaired original strategies conserve budget and capacity", {
  strategies <- c("a_b_n", "a_b_y", "k_b_n", "k_b_y",
                  "a_d_n", "a_d_y", "a_d_a", "k_d_n", "k_d_y", "k_d_a")
  cap <- c(A = 2, B = 4, C = 7, D = 15)
  bg <- c(A = 1, B = 3, C = 9, D = 27)
  for (strategy in strategies) {
    fn <- getFromNamespace(paste0("sample_", strategy), "TaxonSampling")
    for (budget in c(0L, 1L, 3L, 12L, 28L)) {
      set.seed(100 + budget)
      out <- fn(budget, c(A = 0L, B = 0L, C = 0L, D = 0L), cap, bg)
      expect_equal(sum(out), budget, info = strategy)
      expect_true(all(out >= 0 & out <= cap), info = strategy)
      expect_identical(names(out), names(cap), info = strategy)
    }
    expect_error(fn(29L, setNames(integer(4), names(cap)), cap, bg),
                 "capacity", info = strategy)
  }
})

test_that("repaired legacy formulas respect required lower allocations", {
  fn_names <- c("sample_a_b_n", "sample_a_b_y", "sample_k_b_n",
                "sample_k_b_y", "sample_a_d_n", "sample_a_d_y",
                "sample_a_d_a", "sample_k_d_n", "sample_k_d_y", "sample_k_d_a")
  lower <- c(A = 1L, B = 0L, C = 2L)
  cap <- c(A = 2, B = 3, C = 4)
  bg <- c(A = 1, B = 3, C = 5)
  for (name in fn_names) {
    fn <- getFromNamespace(name, "TaxonSampling")
    set.seed(43)
    x <- fn(4L, lower, cap, bg)
    expect_equal(sum(x), 7L, info = name)
    expect_true(all(x >= lower & x <= cap), info = name)
  }
})

test_that("deterministic legacy diversity repeats complete coverage rounds", {
  lower <- c(A = 0L, B = 0L, C = 0L)
  cap <- c(A = 10, B = 10, C = 10)
  bg <- c(A = 1, B = 2, C = 3)
  for (name in c("sample_a_d_n", "sample_k_d_n")) {
    fn <- getFromNamespace(name, "TaxonSampling")
    set.seed(1)
    x <- fn(7L, lower, cap, bg)
    expect_equal(sum(x), 7L)
    expect_equal(sort(unname(x)), c(2L, 2L, 3L))
  }
})

test_that("modern allocator delegates to the repaired legacy policy", {
  fn <- getFromNamespace(".ts_allocate", "TaxonSampling")
  lower <- c(A = 0L, B = 1L, C = 0L)
  cap <- c(A = 2, B = 5, C = 7)
  bg <- c(A = 1, B = 4, C = 9)
  for (method in c("balanced", "diversity")) {
    for (mode in c("no", "yes", "after_first_round")) {
      if (method == "balanced" && mode == "after_first_round") next
      for (sampling in c("agnostic", "known_species")) {
        strategy <- paste0("sample_", if (sampling == "agnostic") "a" else "k",
                           "_", if (method == "balanced") "b" else "d", "_",
                           switch(mode, no = "n", yes = "y", after_first_round = "a"))
        legacy <- getFromNamespace(strategy, "TaxonSampling")
        set.seed(44)
        a <- fn(8L, lower, cap, bg, method, mode, sampling)
        set.seed(44)
        b <- legacy(8L - sum(lower), lower, cap, bg)
        expect_identical(a, b, info = strategy)
      }
    }
  }
})
