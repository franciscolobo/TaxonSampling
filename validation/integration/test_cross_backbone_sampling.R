# Integration smoke test. Run from the TaxonSampling repository root:
# Rscript tests/test_cross_backbone_sampling.R taxonomy_data/taxonomy_db
args <- commandArgs(trailingOnly = TRUE)
db_dir <- if (length(args)) args[[1L]] else 'taxonomy_data/taxonomy_db'
stopifnot(dir.exists(db_dir))
for (pkg in c('DBI', 'duckdb')) stopifnot(requireNamespace(pkg, quietly = TRUE))
source('R/build_cross_identity_db.R')
source('R/get_counts.R')
source('R/run_TS.R')
paths <- setNames(file.path(db_dir, paste0(c('ncbi','ott','synthetic'), '.duckdb')),
                  c('ncbi','ott','synthetic'))
stopifnot(all(file.exists(paths)), file.exists(file.path(db_dir, 'cross_identity.duckdb')))
cross <- file.path(db_dir, 'cross_identity.duckdb')
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = cross, read_only = TRUE)
links <- DBI::dbGetQuery(con, "SELECT source_a, id_a, source_b, id_b
                             FROM identity_links_status WHERE mapping_status='unique'")
transitive <- DBI::dbGetQuery(con, "SELECT ncbi_id, synthetic_id
    FROM transitive_links WHERE mapping_status='unique_path'")
DBI::dbDisconnect(con, shutdown = TRUE)
# Only test a direction if it has at least two distinct uniquely mapped inputs.
pairs <- list(
  list('ncbi','ott',links$id_a[links$source_a=='ncbi' & links$source_b=='ott']),
  list('ott','ncbi',links$id_b[links$source_a=='ncbi' & links$source_b=='ott']),
  list('ott','synthetic',links$id_a[links$source_a=='ott' & links$source_b=='synthetic']),
  list('synthetic','ott',links$id_b[links$source_a=='ott' & links$source_b=='synthetic']),
  list('ncbi','synthetic',transitive$ncbi_id),
  list('synthetic','ncbi',transitive$synthetic_id)
)
for (p in pairs) {
  from <- p[[1]]; to <- p[[2]]
  candidates <- head(unique(as.character(p[[3]])), 100L)
  if (length(candidates) < 2L) stop('Insufficient unique links for ', from, ' -> ', to)
  # Test mapping separately from sampling to avoid relying on candidate IDs
  # being retained in the sampled synthetic tree.
  ids <- data.frame(taxID = candidates, seqID = paste0('seq', seq_along(candidates)))
  x <- get_counts(taxonomy_db = paths[[to]], ids_df = ids,
                  id_type = from, cross_identity_db = cross,
                  missing_ids = 'drop', verbose = FALSE)
  stopifnot(nrow(x$ids_df) >= 2L, all(x$ids_df$node_id %in% x$nodes$id))
  root <- x$metadata$root_id
  result <- run_TS(x, taxon = root, m = 2L, verbose = FALSE)
  stopifnot(length(result$outputIDs) == 2L,
            length(unique(result$outputIDs)) == 2L,
            all(result$outputIDs %in% ids$taxID),
            all(result$outputMapping$node_id %in% x$nodes$id))
  cat('PASS ', from, ' -> ', to, '\n', sep = '')
}
cat('Cross-backbone integration smoke tests passed.\n')
