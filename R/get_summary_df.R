#' Generate a summary dataframe from the result of run_TS
#'
#' Returns a summary data frame containing relevant counts for each taxID
#' associated with the sequences passed to [get_counts()]
#' (parameters `ids_df` or `ids_file`) in the pre-processing stage before the
#' call to [run_TS()].
#'
#' @param taxlist list object of class _taxonsampling_, returned by
#' [run_TS()]
#'
#' @return Data frame containing fields `$taxID`, `$level` (taxonomic level),
#' `$name` (scientific name), `species_count` (total number of species or leaf
#' nodes under that taxonomy ID), `$sample_count` (number of sampled species or
#' leaf nodes under that taxonomy ID) and `$seq_count` (number of available
#' sequences under that taxonomy ID, passed in the call to [get_counts()]).
#'
#' @export
#'

get_summary_df <- function(taxlist){

  assertthat::assert_that(inherits(taxlist, "taxonsampling"))

  # Get ancestry list of (filtered) nodes
  ancestry <- get_taxID_spp_counts(nodes = taxlist$nodes,
                                   start_from_species = FALSE,
                                   what = "ancestry",
                                   verbose = FALSE)


  # Isolate only ancestries of sampled taxIDs
  idx <- which(sapply(ancestry,
                      function(x, ids){x[[1]] %in% ids},
                      ids = taxlist$outputIDs))

  samples <- table(unlist(ancestry[idx]))
  samples <- data.table::data.table(taxID = as.integer(names(samples)),
                                    sample_count = as.integer(samples))

  # Prepare nodes for joining
  nodes <- taxlist$nodes[, -2]
  names(nodes)[1] <- "taxID"

  # Get sequence counts
  seqs <- data.table::data.table(taxID     = as.integer(names(taxlist$countIDs)),
                                 seq_count = as.integer(taxlist$countIDs))

  # Get species counts
  spps <- data.table::as.data.table(taxlist$spp_df)

  # Join all info
  out <- seqs[spps, on = "taxID"]
  out <- samples[out, on = "taxID"]
  out <- nodes[out, on = "taxID"]
  out$sample_count <- ifelse(is.na(out$sample_count), 0, out$sample_count)

  # Filter only ancestries of taxIDs related to the root node
  idx <- which(sapply(ancestry,
                      function(x, tax){(tax %in% x) && (x[[1]] != tax)},
                      tax = taxlist$ts.params$taxon))

  ids <- sapply(ancestry[idx], function(x){x[[1]]})
  out <- out[out$taxID %in% ids, ]

  return(as.data.frame(out[, c(1:3,6:4)]))
}
