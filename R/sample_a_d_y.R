# Corrected original sample_a_d_y allocation strategy.
sample_a_d_y <- function(m, m_i, childrenCount, childrenCountSpp, ...) {
  .ts_legacy_fixed_allocate(m, m_i, childrenCount, childrenCountSpp,
    method = "diversity",
    randomize = "yes",
    sampling = "agnostic")
}
