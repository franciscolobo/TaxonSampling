sample_a_b_y <- function(m, m_i, childrenCount, childrenCountSpp, ...){
  x <- table(sample(names(childrenCount),
                    size = m + sum(m_i),
                    replace = TRUE))

  for (nn in names(m_i)){
    if(!(nn %in% names(x))) {
      x[nn]  <- 0
    }
  }

  while(any(x - m_i < 0)){
    id1 <- which.max(x - m_i)[1]
    id2 <- which.min(x - m_i)[1]
    x[id2] <- x[id2] + 1
    x[id1] <- x[id1] - 1
  }

  return(x)
}
