sample_a_b_y <- function(m, m_i, childrenCount, childrenCountSpp, ...){
  # TODO: why sample with replacement?
  m_i <- table(sample(names(childrenCount), m, replace = TRUE))

  return(m_i)
}
