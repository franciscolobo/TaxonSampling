sample_a_b_n <- function(m, m_i, childrenCount, childrenCountSpp, ...){
  m_i                  <- m_i + floor(m / length(childrenCount))
  sampledChildren      <- sample(names(childrenCount), m - sum(m_i))
  m_i[sampledChildren] <- m_i[sampledChildren] + 1

  return(m_i)
}
