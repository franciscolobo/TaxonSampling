#' Write sequences associated with sampled IDs to a multifasta file.
#'
#' @param taxlist list containing vector `$outputIDs` of sampled IDs (returned
#' by [run_TS()]).
#' @param seq_file character string with the path to the multifasta file
#' containing the input sequences.
#' @param verbose logical: regulates function echoing to console.
#'
#' @return Updated `taxlist` containing field `outputSeqs`, a data frame with
#' information on the sampled sequences.

extract_sequences <- function(taxlist, seq_file, verbose = TRUE) {

  if(verbose) message("Reading sequences from seq file")
  x <- seqinr::read.fasta(file = seq_file, as.string = TRUE)
  x <- x[taxlist$ids_df$seqID[taxlist$ids_df$taxID %in% taxlist$outputIDs]]
  x <- lapply(x, toupper)
  y <- taxlist$ids_df[taxlist$ids_df$seqID %in% names(x), ]
  names(x) <- paste0(y$seqID, "|TaxID:", y$taxID)

  taxlist$outputSeqs <- x
  return(taxlist)
}
