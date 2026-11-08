#' Run Taxon Sampling
#'
#' Run the TaxonSampling method to return a sample of taxonomic IDs according
#' to the desired balance / diversity.
#'
#' @param taxlist list object of class _taxonsampling_, returned by
#' [get_species_count()].
#' @param taxon Taxon ID from which to start sampling children taxa (single
#' character or integer value)
#' @param m desired sample size
#' @param seq_file character string with the path to the multifasta file
#' containing the input sequences.
#' @param out_file character string naming a file to save the output (a
#' multifasta file). Ignored if `seq_file == NULL`.
#' @param method sampling method to use. Accepts "balanced" (favors balanced taxa
#' representation) or "diversity" (favors maximized taxa representation)
#' @param randomize randomization strategy: should the algorithm choose IDs
#' randomly ("yes"), maintaining a balanced allocation ("no"), or with a
#' balanced allocation at the top taxonomic level and randomized afterwards
#' ("after_first_round")?
#' @param replacement logical flag: should the algorithm allow repeated IDs in
#' the output (if needed to reach m IDs in the output with maximized taxonomy
#' diversity).
#' @param ignoreIDs vector (character or integer) of IDs that must not appear in
#' the output.
#' @param requireIDs vector (character or integer) of IDs that must appear
#' in the output. Notice that `ignoreIDs` has precedence over `requireIDs`,
#' i.e., IDs that occur in both will be ignored. `requireIDs` that are children
#' of any `ignoreIDs` will also be ignored.
#' @param sampling sampling mode. Accepts "agnostic" (sample species in a
#' diversity-agnostic manner) or "known_species" (sample based on known species
#' diversity).
#' @param verbose logical: regulates function echoing to console.
#'
#' @return Input object`taxlist` updated with vector `$outputIDs` of
#' sampled IDs and list `$outputSeqs` (if seq_file is not `NULL`) containing
#' information about the sequences sampled.
#'
#' @export

run_TS <- function(taxlist, taxon, m,
                   seq_file         = NULL,
                   out_file         = NULL,
                   method           = "diversity",
                   randomize        = "no",
                   replacement      = FALSE,
                   ignoreIDs        = NULL,
                   requireIDs       = NULL,
                   sampling         = "agnostic",
                   verbose          = TRUE) {

  # ===========================================================================
  # Sanity checks
  assertthat::assert_that(is.list(taxlist),
                          inherits(taxlist, "taxonsampling"),
                          is.character(taxon) || is.numeric(taxon),
                          assertthat::is.count(m),
                          is.null(seq_file) ||
                            (is.character(seq_file) &&
                               length(seq_file) == 1 &&
                               file.exists(seq_file)),
                          is.null(out_file) ||
                            (is.character(out_file) && length(out_file) == 1),
                          is.character(randomize), length(randomize) == 1,
                          is.character(sampling), length(sampling) == 1,
                          is.character(method), length(method) == 1,
                          sampling  %in% c("agnostic", "known_species"),
                          randomize %in% c("yes", "no", "after_first_round"),
                          method    %in% c("diversity", "balanced"),
                          is.logical(replacement), length(replacement) == 1,
                          is.null(ignoreIDs) ||
                            is.numeric(ignoreIDs) ||
                            is.character(ignoreIDs),
                          is.null(requireIDs) ||
                            is.numeric(requireIDs) ||
                            is.character(requireIDs),
                          is.logical(verbose), length(verbose) == 1)

  if (randomize == "after_first_round" && method == "balanced"){
    stop('Combination of randomize == "after_first_round" and method == "balance" not possible')
  }

  # Add all input parameters to taxlist
  taxlist$ts.params         <- as.list(environment())
  taxlist$ts.params$taxlist <- NULL

  # Force all IDs to character
  taxlist$ts.params$taxon      <- as.integer(taxon)
  taxlist$ts.params$ignoreIDs  <- as.integer(ignoreIDs)
  taxlist$ts.params$requireIDs <- as.integer(requireIDs)

  # ===========================================================================

  # Process ignoreIDs and requireIDs
  taxlist <- process_ignoreIDs(taxlist)
  taxlist <- process_requireIDs(taxlist)

  # Reduce node information to the necessary only, reduces search time.
  taxlist$nodes <- taxlist$nodes[taxlist$nodes$id %in% names(taxlist$countIDs), ]

  # Ensure m <= number of valid ids.
  m <- min(m, length(intersect(taxlist$ids_df$taxID, names(taxlist$countIDs))))

  # Call the TS algorithm itself.
  if(verbose) message("\rRunning recursive sampling")
  taxlist$ts.process <- c(taxlist$ts.process,
                          list(taxon = taxon, m = m))
  taxlist$outputIDs  <- ts_recursive(taxlist, verbose)
  if(verbose) cat("\r", paste(rep(" ", 40), collapse = ""), "\r")

  if(!is.null(seq_file)){
    taxlist <- extract_sequences(taxlist, seq_file, verbose)

    if(!is.null(out_file)){
      if(!dir.exists(dirname(out_file))) {
        dir.create(dirname(out_file), recursive = TRUE)
      }
      if(verbose) message("Saving sampled sequences to file: ", out_file)
      seqinr::write.fasta(sequences = taxlist$outputSeqs,
                          names     = names(taxlist$outputSeqs),
                          file.out  = out_file)
    }
  }

  # Remove temporary elements used in processing
  taxlist$ts.process <- NULL

  class(taxlist) <- unique(c(class(taxlist), "taxonsampling"))

  if(verbose) message("Done!")
  return(taxlist)
}
