# Corrected original sample_k_d_y allocation strategy.
sample_k_d_y <- function(m, m_i, childrenCount, childrenCountSpp, ...) {
  .ts_legacy_fixed_allocate(m, m_i, childrenCount, childrenCountSpp,
    method = "diversity",
    randomize = "yes",
    sampling = "known_species")
}
