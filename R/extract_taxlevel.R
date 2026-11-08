#' Extract the desired ancestry level of a list of taxonomy IDs
#'
#' @param taxIDs vector of taxonomy IDs of interest
#' @param taxlevel desired taxonomic level to be extracted (e.g., `"order"`)
#' @param nodes data.frame containing the pre-processed information about
#' the NCBI taxonomy structure. This can be generated, e.g., by using
#' [get_counts()], or as object `$nodes` of a `taxonsampling` object.
#'
#' @return Data frame containing the taxonomy IDs of interest and a column with
#' the corresponding taxonomy IDs of their desired ancestry level.
#'
#' @export
#'
extract_taxlevel <- function(taxIDs, taxlevel, nodes){

  # ===========================================================================
  # Sanity checks
  assertthat::assert_that(is.data.frame(nodes),
                          is.character(taxIDs) || is.integer(taxIDs),
                          is.character(taxlevel),
                          length(taxlevel) == 1)

  taxlevel <- tolower(taxlevel)

  # ===========================================================================

  # Get ancestry list from nodes
  ancestry <- get_taxID_spp_counts(nodes = nodes,
                                   what    ="ancestry",
                                   verbose = FALSE)

  # filter to contain only the ancestry of target taxIDs
  ancestry <- sapply(ancestry,
                     function(x,taxIDs){
                       if(x[[1]] %in% taxIDs) return(x)
                       return(NULL)
                     }, taxIDs = taxIDs)
  ancestry <- ancestry[which(!sapply(ancestry, is.null))]
  names(ancestry) <- as.character(sapply(ancestry, function(x) x[1]))

  # Extract taxlevel of interest and return as data frame
  ancestry <- lapply(ancestry,
                     function(x, nodes, taxlevel){
                       nodes <- data.table::as.data.table(nodes)
                       x <- data.table::data.table(id = x)
                       x <- nodes[x, on = "id"][, c(1,3)]
                       x <- x[tolower(x$level) == tolower(taxlevel), ]
                     }, nodes = nodes, taxlevel = taxlevel)

  ancestry <- data.table::rbindlist(ancestry, use.names = TRUE, idcol = TRUE)
  names(ancestry)[1:2] <- c("taxID", taxlevel)

  return(as.data.frame(ancestry[, -3]))

}





