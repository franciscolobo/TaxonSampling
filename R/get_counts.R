#' Return the relevant counts for running TaxonSampling
#'
#' This function evaluates the taxonomical structure from a set of input IDs
#' and returns the count of each taxon ID with non-zero occurrences in the
#' taxonomy graph, as well as the count of known species under each taxon.
#'
#' @param taxonomy_path path to folder containing the NCBI taxonomy files
#' (i.e., the extracted contents of _taxdump.zip_, which can be downloaded from
#' <ftp://ftp.ncbi.nih.gov/pub/taxonomy/> or retrieved using
#' [retrieve_NCBI_taxonomy()]).
#' @param ids_file path to a tab-separated file containng two two columns,
#' with the input taxon IDs in the first
#' column, and the corresponding sequence IDs (corresponding
#' to the ID strings int the multifasta input file, without
#' the ">" line starter) in the second column.
#' @param ids_df two-column data frame with the input taxon IDs in the first
#' column, and the corresponding sequence IDs (corresponding
#' to the ID strings int the multifasta input file, without
#' the ">" line starter) in the second column. Ignored if
#' `ids_file` is not `NULL`.
#' @param spp_df two-column data frame with the input taxon IDs in the first
#' column, and the corresponding number of known species in the second column.
#' This can be generated using [get_taxID_spp_counts()].
#' Ignored if `spp_file` is not `NULL`.
#' @param spp_file path to a tab-separated file containing two two columns,
#' with the input taxon IDs in the first column, and the corresponding number of
#' known species in the second column. If both `spp_file` and `spp_df` are `NULL`, then
#' the counts are computed internally based either on the
#' taxonomy files under `taxonomy_path`.
#' @param verbose logical: regulates function echoing to console.
#' @param start_from_species logical, passed down to [get_taxID_spp_counts()]
#' (only if `spp_file` and `spp_file` are `NULL`).
#'
#' @return list object of class _taxonsampling_, containing:
#' \itemize{
#'     \item `$ids_df`: data.frame with taxon IDs in column 1 and corresponding
#'     sequence IDs in column 2, as loaded from `ids_file` or passed directly as
#'     input. Filtered to maintain only IDs that exist in `$nodes` and to remove
#'     duplicated IDs.
#'     \item `$nodes`: data.frame containing the pre-processed information about
#'     the NCBI taxonomy structure, extracted from file _nodes.dmp_ of the
#'     taxonomy files. Filtered to keep only nodes with IDs present in `$ids_df`
#'     and with a non-zero total count.
#'     \item `$countIDs`: numeric vector with the counts of the number of
#'     taxonomy nodes (of all levels) under each taxon ID.
#'     \item `spp_df`: two-column data frame with the input taxon IDs in the
#'     first column, and the corresponding number of known species in the second
#'    column. Filtered to have only the IDs present in `$countIDs`.
#' }
#'
#' @export
#'
#' @importFrom data.table :=

get_counts <- function(taxonomy_path,
                       ids_file      = NULL,
                       ids_df        = NULL,
                       spp_file      = NULL,
                       spp_df        = NULL,
                       start_from_species = FALSE,
                       verbose       = TRUE) {

  # ===========================================================================
  # Sanity checks
  assertthat::assert_that(is.null(ids_df) ||
                            (is.data.frame(ids_df) && nrow(ids_df) > 0),
                          is.null(ids_file) ||
                            (is.character(ids_file) && length(ids_file) == 1),
                          (is.null(ids_df) + is.null(ids_file)) < 2,
                          is.character(taxonomy_path),                               length(taxonomy_path) == 1 &&
                            dir.exists(taxonomy_path),
                          is.null(spp_df) || (
                            is.data.frame(spp_df) &&
                              ncol(spp_df) >= 2 &&
                              nrow(spp_df) > 0),
                          is.null(spp_file) ||
                            (is.character(spp_file) &&
                               length(spp_file) == 1 &&
                               file.exists(spp_file)),
                          is.logical(verbose),
                          length(verbose) == 1,
                          is.logical(start_from_species),
                          length(start_from_species) == 1)

  # ===========================================================================
  # Load required files

  odt <- options("datatable.showProgress")
  options(datatable.showProgress = FALSE)
  # Load ids from file if required
  if(!is.null(ids_file)) {
    if(file.exists(ids_file)){
      if(verbose) message('Reading ids file')
      ids_df <- as.data.frame(
        data.table::fread(ids_file, sep = "\t",
                          colClasses = c("integer", "character"),
                          col.names = c("taxID", "seqID"),
                          verbose = FALSE))
    } else {
      stop("File ", ids_file, " not found.")
    }
  } else {
    names(ids_df) <- c("taxID", "seqID")
    ids_df$taxID <- as.integer(ids_df$taxID)
  }

  ids_df$seqID <- gsub("\\t", "", ids_df$seqID)


  # Load nodes from file
  if(verbose) message('Reading nodes file')
  nodes <- data.table::fread(paste(taxonomy_path, "nodes.dmp", sep = "/"),
                             sep = "|", strip.white = TRUE,
                             colClasses = c("integer", "integer", "character", rep("NULL", 16)),
                             col.names = c("id", "parent", "level"),
                             verbose = FALSE)
  nodes$level <- gsub("\\t", "", nodes$level)

  if(verbose) message('Reading names file')
  name <- NULL
  names <- data.table::fread(paste(taxonomy_path, "names.dmp", sep = "/"),
                             sep = "|", strip.white = TRUE,
                             colClasses = c("integer", "character", "NULL", "character", "NULL"),
                             col.names = c("id", "name", "status"),
                             verbose = FALSE)
  names[, 2:3] <- lapply(names[, 2:3], function(x) gsub("\\t", "", x))
  names <- names[names$status == "scientific name", -c("status")]
  nodes[names, name := name, on = c("id")]


  # Load ids from file if required
  if(!is.null(spp_file)){
    if(verbose) message('Reading spp file')
    spp_df <- as.data.frame(
      data.table::fread(spp_file, sep = "\t",
                        colClasses = c("integer", "integer"),
                        col.names = c("taxID", "species_count"),
                        verbose = FALSE))
  } else if (!is.null(spp_df)){
    names(spp_df)[1:2] <- c("taxID", "species_count")
    spp_df$taxID <- as.integer(spp_df$taxID)
    spp_df$species_count <- as.integer(spp_df$species_count)
  } else {
    spp_df <- get_taxID_spp_counts(taxonomy_path,
                                   nodes   = nodes,
                                   verbose = verbose,
                                   start_from_species = start_from_species)
  }

  # Return data.table options to previous state
  options(odt)

  # ===========================================================================
  # Process taxonomy counts

  # Filter IDs that aren't part of NCBI notation.
  idx <- which(!(ids_df$taxID %in% nodes$id))
  if (length(idx) > 0) {
    warning("The following IDs are not found in the NCBI taxonomy files and will be ignored:\n",
            paste(ids_df$taxID[idx], collapse = "\n"))
    ids_df <- ids_df[-idx, ]
  }

  # Remove duplicates
  idx <- which(duplicated(ids_df$taxID))
  if (length(idx) > 0) {
    warning("Some IDs are duplicated, only the first occurrence will be used.")
    ids_df <- ids_df[-idx, ]
  }

  # Count how often each taxonomy ID occurs in the dataset
  countIDs            <- rep.int(0, nrow(nodes))
  names(countIDs)     <- nodes$id
  nodes_parent        <- nodes$parent
  names(nodes_parent) <- nodes$id
  searchIDs           <- ids_df$taxID
  countIDs[as.character(searchIDs)] <- 1

  if(verbose) message("Counting taxonomy IDs")
  while (length(searchIDs) > 0) {
    searchIDs <- nodes_parent[as.character(searchIDs)]
    parentage <- table(searchIDs)
    countIDs[names(parentage)] <- countIDs[names(parentage)] + parentage
    searchIDs <- searchIDs[searchIDs != 1]
    if(verbose) cat("\r", sprintf("--> %06d taxIDs still being counted.", length(searchIDs)))
  }
  if(verbose) cat("\r", paste(rep(" ", 40), collapse = ""))

  countIDs <- countIDs[countIDs > 0]

  # ===========================================================================
  # Filter only the species counts of taxa listed in countIDs
  spp_df   <- spp_df[spp_df$taxID %in% names(countIDs), ]

  # Simplify nodes to reduce overhead
  nodes <- nodes[nodes$id %in% as.numeric(names(countIDs)), ]

  taxlist <- list(nodes    = nodes,
                  ids_df   = ids_df,
                  countIDs = countIDs,
                  spp_df   = spp_df)

  class(taxlist) <- c(class(taxlist), "taxonsampling")

  return(taxlist)
}
