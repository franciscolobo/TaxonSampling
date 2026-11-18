# In-memory fixture compatible with the current run_TS() interface.
# Nodes may use arbitrary, nonnumeric native backbone identifiers.
make_ts_fixture <- function(prefix = "ott:", n_per_tip = 3L,
                            direct_root = FALSE, bg = NULL,
                            count_mode = "input_ids") {
  stopifnot(n_per_tip >= 1L)
  keys <- c("root", "A", "B", "A1", "A2", "B1", "B2")
  ids <- paste0(prefix, keys)
  parent <- paste0(prefix, c("root", "root", "root", "A", "A", "B", "B"))
  nodes <- data.frame(id = ids, parent = parent, name = keys,
                      rank = c("no rank", "phylum", "phylum", rep("species", 4L)),
                      depth = c(0L, 1L, 1L, rep(2L, 4L)),
                      child_count = c(2L, 2L, 2L, rep(0L, 4L)),
                      stringsAsFactors = FALSE)
  leaves <- paste0(prefix, c("A1", "A2", "B1", "B2"))
  obs <- do.call(rbind, lapply(seq_along(leaves), function(i) {
    data.frame(taxID = paste0("tax_", i, "_", seq_len(n_per_tip)),
               seqID = paste0("seq_", i, "_", seq_len(n_per_tip)),
               node_id = leaves[i], stringsAsFactors = FALSE)
  }))
  if (direct_root) obs <- rbind(obs, data.frame(taxID = "tax_root",
    seqID = "seq_root", node_id = ids[1L], stringsAsFactors = FALSE))
  if (is.null(bg)) bg <- c(12, 6, 6, 3, 3, 3, 3)
  stopifnot(length(bg) == nrow(nodes))
  result <- list(nodes = nodes, ids_df = obs,
                 countIDs = data.frame(taxID = ids, count = 0L),
                 spp_df = data.frame(taxID = ids, species_count = as.numeric(bg)),
                 metadata = list(count_mode = count_mode, root_id = ids[1L]))
  class(result) <- "taxonsampling"
  result
}

sample_fixture <- function(x = make_ts_fixture(), m = 4L, ...) {
  TaxonSampling::run_TS(x, taxon = x$metadata$root_id,
                        m = m, verbose = FALSE, ...)
}

# Node ancestry and clade coverage, independent of the algorithm under test.
fixture_descendants <- function(nodes, ancestor) {
  result <- ancestor
  repeat {
    new <- nodes$id[nodes$parent %in% result & !nodes$id %in% result]
    if (!length(new)) break
    result <- c(result, new)
  }
  result
}

covered <- function(result, ancestor) {
  any(result$outputMapping$node_id %in%
        fixture_descendants(result$nodes, ancestor))
}
