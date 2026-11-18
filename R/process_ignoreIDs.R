process_ignoreIDs <- function(taxlist) {
  ignoreIDs <- taxlist$ts.params$ignoreIDs
  if (!is.null(ignoreIDs) & length(ignoreIDs) > 0) {
    # Filter IDs that aren't part of NCBI notation.
    idx <- which(!(ignoreIDs %in% taxlist$nodes$id))
    if (length(idx) > 0) {
      warning("The following IDs are not found in the NCBI taxonomy files and will be ignored:\n",
              paste(ignoreIDs[idx], collapse = "\n"))
      ignoreIDs <- ignoreIDs[-idx, ]
    }

    # Remove duplicates
    idx <- which(duplicated(ignoreIDs))
    if (length(idx) > 0) {
      warning("Some IDs are duplicated, only the first occurrence will be used.")
      ignoreIDs <- ignoreIDs[-idx]
    }

    # Make sure none of the ignoreIDs is a parent/ancestor of another.
    # It simplifies the processing of the following steps.
    for (id in ignoreIDs) {
      if (any(ignoreIDs[ignoreIDs != id] %in% CHNOSZ::allparents(id, nodes = taxlist$nodes))) {
        ignoreIDs <- ignoreIDs[ignoreIDs != id]
      }
    }

    # Subtract the ignoreIDs from the countIDs.
    for (id in ignoreIDs) {
      parentIDs <- as.character(CHNOSZ::allparents(id, nodes = taxlist$nodes))
      taxlist$countIDs[parentIDs] <- taxlist$countIDs[parentIDs] - taxlist$countIDs[id]
    }
    taxlist$countIDs <- taxlist$countIDs[taxlist$countIDs > 0]

    # Prune out the children of removed nodes.
    orphans <- taxlist$nodes$id[!(taxlist$nodes$parent %in% names(taxlist$countIDs)) &
                                  (taxlist$nodes$id %in% names(taxlist$countIDs))]
    while (length(orphans) > 0) {
      taxlist$countIDs <- taxlist$countIDs[setdiff(names(taxlist$countIDs), orphans)]
      orphans <- taxlist$nodes$id[!(taxlist$nodes$parent %in% names(taxlist$countIDs)) &
                            (taxlist$nodes$id %in% names(taxlist$countIDs))]
    }

  }

  return(taxlist)
}
