# Corrected original sample_a_b_n allocation strategy.
sample_a_b_n <- function(m, m_i, childrenCount, childrenCountSpp, ...) {
  .ts_legacy_fixed_allocate(m, m_i, childrenCount, childrenCountSpp,
    method = "balanced",
    randomize = "no",
    sampling = "agnostic")
}
