#' Extract exactly the sampled sequences, in sampling order.
#' @export
extract_sequences <- function(taxlist, seq_file, verbose = TRUE) {
  if (is.null(taxlist$outputMapping))
    stop("Run run_TS() before extract_sequences().")
  if (!requireNamespace("seqinr", quietly = TRUE))
    stop("Package seqinr is required for FASTA extraction.")
  if (verbose) message("Reading FASTA: ", seq_file)
  fasta <- seqinr::read.fasta(file = seq_file, as.string = TRUE,
                              seqtype = "DNA", whole.header = FALSE)
  requested <- as.character(taxlist$outputMapping$seqID)
  missing <- unique(requested[!requested %in% names(fasta)])
  if (length(missing)) stop("FASTA records missing: ", paste(missing, collapse = ", "))
  sequences <- lapply(fasta[requested], function(x) toupper(as.character(x)))
  names(sequences) <- make.unique(paste0(requested, "|TaxID:",
                                             taxlist$outputMapping$taxID), sep = "_repeat")
  taxlist$outputSeqs <- sequences
  taxlist
}
