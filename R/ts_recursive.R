# Receives a group of Taxonomy IDs and the desired sample size m.
# Returns a sampled vector with maximized taxonomy diversity.
# Assumes that every input id is unique.
# May return repeated IDs if replacement == TRUE, to increase
# diversity when method == "diversity".
#
# Args:
# list containing fields:
#   taxlist$ts.params$taxon:  (char/integer) Taxon from which to start sampling
#                             children taxa.
#   taxlist$ts.params$m:      (integer) desired sample size.
#   taxlist$nodes:            (data.frame) pre-processed information from NCBI
#                             taxonomy.
#   taxlist$ts.params$countIDs: (numeric) count of taxnomoy IDs belonging to
#                             each taxon, created by get_taxonomy_counts()
#   replacement:              (logical) whether the algorithm allows to repeat
#                             IDs in order to maximize taxonomy diversity and to
#                             reach m IDs in the output.
#   randomize:                (char) whether the algorithm will choose IDs
#                             randomly or maintaining a balanced allocation (m_i
#                             differing by at most 1 if the maximum possible
#                             value wasn't reached).
#   method:                   (char) whether it favors balanced taxa
#                             representation ("balance") or maximized taxa
#                             representation ("diversity").
#   requireIDs:               (char) IDs that must appear in the output.
#
# Returns:
#   updated taxlist containing field:
#   outputIDs:                (vector) vector of IDs with maximized taxonomy
#                             balance or diversity.

ts_recursive <- function(taxlist, verbose = TRUE) {


  # extract relevant variables for recursive call
  taxon <- taxlist$ts.process$taxon
  m     <- taxlist$ts.process$m

  countSpp <- taxlist$spp_df$species_count
  names(countSpp) <- taxlist$spp_df$taxID

  if(verbose) {
    # This isn't intended to mean much, it just gives the user some
    # visual feedback
    cat("\r", paste(rep(" ", 80), collapse = ""),
        "\r--> Recursive sampling:",
        paste(rep(".",
                  min(30, ceiling(60 * m/taxlist$ts.params$m))),
              collapse = ""))
  }

  # Sanity check
  if (m <= 0) stop("m less or equal than zero during recursion.")

  # First step: Find the sub-taxa (children nodes) of the current taxon
  # that has members in user-provided data
  children <- taxlist$nodes$id[taxlist$nodes$parent == taxon & taxlist$nodes$id != taxon]
  children <- intersect(children, names(taxlist$countIDs))

  # Condition to end recursion
  if (length(children) == 0) {
    if(taxlist$ts.params$replacement){
      return(rep(taxon, m))
    } else {
      return(taxon)
    }
  }

  childrenCount    <- taxlist$countIDs[as.character(children)]
  childrenCountSpp <- countSpp[as.character(children)]

  # Sanity check
  # In a few cases, the number of sequences may be greater than the number of
  # species. We saw this happening only in terminal nodes where there's both a
  # species and a subspecies, as our mapping of taxID2species ends at the
  # species level.
  # If that's the case, replace the number of known species by the number of
  # child sequences.
  childrenCountSpp[childrenCount > childrenCountSpp] <- childrenCount[childrenCount > childrenCountSpp]


  # For cases when one taxon isn't a leaf node, but is an input ID that should
  # be available to be sampled along with its children.
  # Example of this case is an input with Homo sapiens (9606),
  # H. sapiens neanderthalensis (63221) and H. sapiens ssp. denisova (741158).
  if (sum(childrenCount) < taxlist$countIDs[as.character(taxon)]) {
    childrenCount    <- c(childrenCount, taxlist$countIDs[as.character(taxon)])
    childrenCountSpp <- c(childrenCountSpp, taxlist$countIDs[as.character(taxon)])
    childrenCountSpp[as.character(taxon)] <- 1
    childrenCount[as.character(taxon)]    <- 1
  }

  # Allocation function to use:
  fname <- paste("sample",
                 substr(taxlist$ts.params$sampling, 1, 1),
                 substr(taxlist$ts.params$method, 1, 1),
                 substr(taxlist$ts.params$randomize, 1, 1),
                 sep = "_")

  # Ensure requireIDs are present in the sample in the required quantities
  m_i  <- 0 * childrenCount
  reqs <- which(names(m_i) %in% names(taxlist$ts.process$requireIDs))
  if(length(reqs) > 0){
    m_i[reqs] <- taxlist$ts.process$requireIDs[names(m_i)[reqs]]
  }
  m <- m - sum(m_i)

  # Perform sampling at current recursion level
  m_i <- do.call(fname,
                 args = list(m                = m,
                             m_i              = m_i,
                             childrenCount    = childrenCount,
                             childrenCountSpp = childrenCountSpp))

  outputIDs <- integer()
  for (id in names(m_i)) {
    if (m_i[id] > 0){
      if (id == taxon){
        outputIDs <- c(outputIDs, as.integer(id))
      } else {
        taxlist$ts.process$m     <- m_i[id]
        taxlist$ts.process$taxon <- as.integer(id)
        outputIDs <- c(outputIDs, ts_recursive(taxlist, verbose))
      }
    }
  }

  return(outputIDs)
}
