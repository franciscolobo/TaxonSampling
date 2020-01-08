#-----------------------------------------------------------------------------#
# Parameters for the analysis, uncomment and edit them if you prefer to leave
# it automated for future uses.


#where should we start looking?
root_taxon <- 1


#taxon ids to sequence names
#idsFile <- "inputs/test1"
idsFile <- "data/validation/taxonIDs_2_sequenceIDs.txt"

#fasta files
#multifasta <- "inputs/multifasta"
multifasta <- "data/validation/sequences_with_taxonIDs.fasta"

#NCBI taxonomy files
taxondir <- "taxdump/"

#number of sequences to sample
m <- 200

ignoreIDs <- NULL

requireIDs <- NULL

#requireIDs <- c(2026169, 8364, 57393, 241292, 61967)

ignoreNonLeafID <- NULL

outFile <- "output.fasta"

#-----------------------------------------------------------------------------#

source("bin/TaxonSampling/TaxonSampling.R")

#library("CHNOSZ")
#library("ape")

nodes <- suppressMessages(getnodes(taxondir))

countIDs <- TS_TaxonomyData(idsFile, nodes)

##These variables are highly experimental and probably
#will not to work, hope future versions of TS will make
#them work

#method to use (either 'diversity' or 'balance')
method <- "diversity"

#wheter to randomize
randomize <- "yes"

#allow TS to repeat IDs if needed? ('no' is better to get a higher diversity)
replacement <- "yes"

nodes <- Simplify_Nodes(nodes, countIDs)

outputIDs1 <- list()

if (!is.null(requireIDs)) {
  total_requiredIDs <- length(requireIDs)
  if (!all(is.element(requireIDs, nodes[,1]))) {
    stop("One of the required IDs provided is not present in database")
  }
  outputIDs1 <- intersect(requireIDs, nodes[,1])
  m <- m - total_requiredIDs
}

requireIDs <- NULL

outputIDs2 <- TS_Algorithm(root_taxon, m, nodes, countIDs, method, randomize,
                          replacement, ignoreIDs, requireIDs,
                          ignoreNonLeafID)

outputIDs <- c(outputIDs1, outputIDs2)

#print(requireIDs)

WriteFasta(idsFile, multifasta, outputIDs, outFile)
