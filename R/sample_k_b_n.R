sample_k_b_n <- function(m, m_i, childrenCount, childrenCountSpp, ...){
  m_i                  <- floor(m / length(childrenCount))
  sampledChildren      <- sample(names(childrenCount),
                                 size = m - sum(m_i),
                                 prob = childrenCountSpp / sum(childrenCountSpp))
  m_i[sampledChildren] <- m_i[sampledChildren] + 1

  return(m_i)
}