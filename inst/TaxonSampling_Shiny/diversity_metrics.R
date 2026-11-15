# Diversity and representativeness metrics for mutually exclusive rank-level groups.
# No external dependencies; source from app.R.
ts_hill <- function(p, q) {
  p <- p[is.finite(p) & p > 0]
  if (!length(p)) return(NA_real_)
  p <- p / sum(p)
  if (q == 0) return(as.numeric(length(p)))
  if (q == 1) return(exp(-sum(p * log(p))))
  if (q == 2) return(1 / sum(p^2))
  (sum(p^q))^(1/(1-q))
}

ts_jsd <- function(p, q) {
  if (length(p) != length(q) || any(!is.finite(p)) ||
      any(!is.finite(q)) || any(p < 0) || any(q < 0) ||
      sum(p) <= 0 || sum(q) <= 0) return(NA_real_)
  p <- p / sum(p); q <- q / sum(q); m <- (p + q)/2
  kl_half <- function(x) sum(x[x > 0] * log2(x[x > 0] / m[x > 0]))
  max(0, min(1, (kl_half(p) + kl_half(q))/2))
}

ts_diversity_metrics <- function(d, rank_label) {
  stopifnot(all(c('background','input','sampled') %in% names(d)))
  counts <- as.matrix(d[, c('background','input','sampled'), drop=FALSE])
  storage.mode(counts) <- 'double'
  if (any(!is.finite(counts)) || any(counts < 0)) stop('Invalid taxonomic counts.')
  totals <- colSums(counts)
  if (any(totals <= 0)) stop('Reference, input and sampled datasets must have positive counts at the selected rank.')
  freq <- sweep(counts, 2, totals, '/')
  colnames(freq) <- c('Reference','Biased database','Sampled set')
  ref <- freq[,1]
  jsd_db <- ts_jsd(ref, freq[,2]); jsd_sample <- ts_jsd(ref, freq[,3])
  reduction <- if (is.finite(jsd_db) && jsd_db > 0)
    100 * (1 - jsd_sample/jsd_db) else NA_real_
  data.frame(
    rank=rank_label, dataset=colnames(freq),
    groups_present=colSums(counts > 0),
    hill_q0=apply(freq, 2, ts_hill, q=0),
    hill_q1=apply(freq, 2, ts_hill, q=1),
    hill_q2=apply(freq, 2, ts_hill, q=2),
    jsd_vs_reference=c(0,jsd_db,jsd_sample),
    reference_group_coverage=colSums(counts > 0 & counts[,1] > 0) / sum(counts[,1] > 0),
    jsd_reduction_pct=c(NA_real_, NA_real_, reduction),
    stringsAsFactors=FALSE
  )
}
