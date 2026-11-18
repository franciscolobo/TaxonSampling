# Corrected original sample_k_d_n allocation strategy.
sample_k_d_n <- function(m, m_i, childrenCount, childrenCountSpp, ...) {
  .ts_legacy_fixed_allocate(m, m_i, childrenCount, childrenCountSpp,
    method = "diversity",
    randomize = "no",
    sampling = "known_species")
}
