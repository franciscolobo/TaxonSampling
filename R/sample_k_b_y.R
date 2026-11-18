# Corrected original sample_k_b_y allocation strategy.
sample_k_b_y <- function(m, m_i, childrenCount, childrenCountSpp, ...) {
  .ts_legacy_fixed_allocate(m, m_i, childrenCount, childrenCountSpp,
    method = "balanced",
    randomize = "yes",
    sampling = "known_species")
}
