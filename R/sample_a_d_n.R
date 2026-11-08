sample_a_d_n <- function(m, m_i, childrenCount, childrenCountSpp, ...){

  while (m > 0 & sum(childrenCount > m_i) <= m) {
    child <- names(childrenCount)[childrenCount > m_i]
    m_i[child] <- m_i[child] + 1
    m <- m - length(child)
  }
  if(sum(m_i) < m){
    child <- sample(names(childrenCount)[childrenCount > m_i], m)
    m_i[child] <- m_i[child] + 1
  }

  return(m_i)
}
