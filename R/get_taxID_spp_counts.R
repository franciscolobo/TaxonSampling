#' Count number of species under each taxonomy ID
#'
#' Extract the species count for each taxonomy ID in the _nodes_ structure
#' of the NCBI taxonomy.
#'
#' @param taxonomy_path path to folder containing the NCBI taxonomy files
#' (i.e., the extracted contents of _taxdump.zip_, which can be downloaded from
#' <ftp://ftp.ncbi.nih.gov/pub/taxonomy/> or retrieved using
#' [retrieve_NCBI_taxonomy()]). Ignored if `nodes` is not `NULL`.
#' @param nodes data.frame containing the pre-processed information about
#' the NCBI taxonomy structure. This can be generated, e.g., by using
#' [CHNOSZ::getnodes()].
#' @param out_file path to file for saving the output
#' (saved as a tsv file). Ignored if `NULL`.
#' @param start_from_species logical. If `TRUE` the counting starts at
#' the 'species' level - i.e., taxon IDs below species (e.g., 'subspecies')
#' receive a count of zero and are not included in the other counts. If
#' `FALSE`, counting starts at the terminal nodes of the taxonomy, regardless of
#' taxonomic level).
#' @param what what to calculate? Accepts `"spp_counts"` for returning the
#' count of species / leaf nodes per taxon ID, or `"ancestry"` for returning
#' the ancestry list.
#' @param verbose logical: regulates function echoing to console.
#'
#' @return either a data frame with species counts per taxon ID
#' or a list with the ancestry of each taxonomic ID.
#'
#' @importFrom data.table := .N
#'
#' @export

get_taxID_spp_counts <- function(taxonomy_path = NULL,
                                 nodes    = NULL,
                                 out_file = NULL,
                                 start_from_species = FALSE,
                                 what    = "spp_counts",
                                 verbose = TRUE){

  # ===========================================================================
  # Sanity checks
  assertthat::assert_that(is.null(taxonomy_path) ||
                            (is.character(taxonomy_path) &&
                               length(taxonomy_path) == 1 &&
                               dir.exists(taxonomy_path)),
                          is.null(nodes) || is.data.frame(nodes),
                          is.null(nodes) + is.null(taxonomy_path) < 2,
                          is.null(out_file) ||
                            (is.character(out_file) &&
                               length(out_file) == 1),
                          is.logical(start_from_species),
                          length(start_from_species) == 1,
                          is.character(what), length(what) == 1,
                          what %in% c("spp_counts", "ancestry"),
                          is.logical(verbose), length(verbose) == 1)

  # ===========================================================================
  # Get nodes table
  if(is.null(nodes)){
    if(verbose) message('Reading nodes file')
    # read taxonomy file
    odt <- options("datatable.showProgress")
    options(datatable.showProgress = FALSE)
    nodes <- data.table::fread(paste(taxonomy_path, "nodes.dmp", sep = "/"),
                               sep = "|", strip.white = TRUE,
                               colClasses = c("integer", "integer", "character", rep("NULL", 16)),
                               col.names = c("id", "parent", "level"),
                               verbose = FALSE)
    options(odt)
  } else {
    nodes <- data.table::as.data.table(nodes[, 1:3])
    names(nodes) <- c("id", "parent", "level")
    nodes$id <- as.integer(nodes$id)
    nodes$parent <- as.integer(nodes$parent)
  }

  nodes$level <- gsub("\\t", "", nodes$level)

  # Add fictional super-root node "0" (makes it easier to remove later)
  nodes$parent[nodes$id == 1] <- 0L
  nodes <- rbind(nodes[1, ], nodes)
  nodes[1, 1:2] <- 0L

  # Initialise data.table variable names to prevent warnings
  id     <- NULL
  parent <- NULL

  if(verbose) message('Extracting ancestry structure')

  # initialise ancestry table
  ids <- nodes[-1, 1:2]
  if(start_from_species) ids <- ids[nodes$level[-1] == 'species', ]

  while(any(ids$parent != 0L)){
    if(verbose) cat(sprintf('\r--> Remaining IDs: %10d', sum(ids$parent != 0L)))
    # Rename columns to facilitate processing (our query is always called 'id')
    names(ids) <- c(paste0("o", (ncol(ids)-1):1), "id")
    # join the parents of our (updated) 'id' column
    ids[nodes, on = list(id), parent := list(parent)]
  }
  ids$parent <- NULL

  if(verbose) cat("\r", paste(rep(" ", 50), collapse = ""), "\r")

  if (what == "spp_counts"){
    # Count occurrences of each taxonID in the ancestry table
    if(verbose) message('Extracting species counts')
    X    <- unname(unlist(ids))
    X    <- data.table::data.table(x = X[X != 0L])
    out <- X[, .N, by = c("x")]
    names(out) <- c("taxID", "species_count")
    out$taxID  <- as.integer(out$taxID)

    # If start_from_species is TRUE, add remaining IDs with a count of zero.
    toadd <- nodes$id[!(nodes$id %in% c(0L, out$taxID))]
    out <- as.data.frame(
      rbind(out,
            data.frame(taxID = toadd,
                       species_count = integer(length(toadd)))))
  } else {
    if(verbose) message('Assembling ancestry list')
    ids[ids == 0L] <- NA
    out <- data.table::transpose(as.list(ids))
    out <- lapply(out, function(x) x[!is.na(x)])
    # names(out) <- as.character(ids[, 1])
  }


  # Save to file
  if(!is.null(out_file)){
    if(!dir.exists(dirname(out_file))) {
      dir.create(dirname(out_file), recursive = TRUE)
    }
    utils::write.table(out, file = out_file, sep = "\t", row.names = FALSE)
  }

  return(out)

}
