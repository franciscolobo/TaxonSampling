# Corrected original sample_a_d_n allocation strategy.
sample_a_d_n <- function(m, m_i, childrenCount, childrenCountSpp, ...) {
  .ts_legacy_fixed_allocate(m, m_i, childrenCount, childrenCountSpp,
    method = "diversity",
    randomize = "no",
    sampling = "agnostic")
}
