process_requireIDs <- function(taxlist) {

  taxlist$ts.process$requireIDs <- numeric()
  requireIDs <- taxlist$ts.params$requireIDs
  if (!is.null(requireIDs) && length(requireIDs) > 0) {
    # Filter IDs that aren't part of NCBI notation.
    idx <- which(!(requireIDs %in% taxlist$nodes$id))
    if (length(idx) > 0) {
      warning("The following IDs are not found in the NCBI taxonomy files and will be ignored:\n",
              paste(requireIDs[idx], collapse = "\n"))
      requireIDs <- requireIDs[-idx, ]
    }

    # Remove duplicates
    idx <- which(duplicated(requireIDs))
    if (length(idx) > 0) {
      warning("Some IDs are duplicated, only the first occurrence will be used.")
      requireIDs <- requireIDs[-idx]
    }

    # Remove requireIDs that are parents of other requireIDs, since this is
    # redundant - EXCEPT if the parent id is at level "species" or below.
    for (id in requireIDs) {
      idpars <- CHNOSZ::allparents(id, nodes = taxlist$nodes)
      idpars <- taxlist$nodes[taxlist$nodes$id %in% idpars, ]
      idpars <- idpars$id[!(idpars$level %in% c("species",
                                                "subspecies",
                                                "varietas"))]
      idpars <- idpars[idpars != id]
      if (length(idpars) > 0) {
        requireIDs <- requireIDs[!(requireIDs %in% idpars)]
      }
    }

    if (length(requireIDs) > 0) {
      parentIDs <- table(unlist(lapply(requireIDs,
                                       CHNOSZ::allparents,
                                       nodes = taxlist$nodes)))
      taxlist$ts.process$requireIDs <- as.numeric(parentIDs)
      names(taxlist$ts.process$requireIDs) <- names(parentIDs)

    }
  }

  return(taxlist)
}
