# Corrected original sample_a_d_a allocation strategy.
sample_a_d_a <- function(m, m_i, childrenCount, childrenCountSpp, ...) {
  .ts_legacy_fixed_allocate(m, m_i, childrenCount, childrenCountSpp,
    method = "diversity",
    randomize = "after_first_round",
    sampling = "agnostic")
}
