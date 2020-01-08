# For parallel processing. for a serial run, do "cores <- 1"
suppressMessages(library("foreach"))
suppressMessages(library("doParallel"))
library("ggplot2")
source("bin/TaxonSampling/TaxonSampling.R")

#Number of bootstraps
n <- 100

#x axis values when plotting results
x <- c(50, 100, 150, 200, 250, 300, 350, 400)

#where should we start looking?
root_taxon <- 1

#taxon ids to sequence names
idsFile <- "data/validation/taxonIDs_2_sequenceIDs.txt"

#fasta files
multifasta <- "data/validation/sequences_with_taxonIDs.fasta"

#NCBI taxonomy files
taxondir <- "taxdump/"

#loading node structure from NCBI Taxonomy
nodes <- suppressMessages(getnodes(taxondir))

#number of taxonomic IDs per node
countIDs <- TS_TaxonomyData(idsFile, nodes)


cores <- 3
if (cores > 1) {
  cl <- makeCluster(cores)
  registerDoParallel(cl)
}
#
comb <- function(x, ...) {
  lapply(seq_along(x),
    function(i) c(x[[i]], lapply(list(...), function(y) y[[i]])))
}

oper <- foreach(i=1:10, .combine='comb', .multicombine=TRUE,
                .init=list(list(), list())) %dopar% {
  list(i+2, i+3)
}

foreach (i = 1:n, .export = c("TS_Algorithm", "RandomSampling",
                              "Evaluate_TS", "TS_TaxonomyData",
                              "TS_Algorithm")) %dopar% {
                              }

countIDs <- TS_TaxonomyData(idsFile, nodes)

# Reduce the node information to the necessary only, reduces search time.
nodes <- nodes[is.element(nodes$id, names(countIDs)), 1:2]

randomize <- "no"

method <- "diversity"
outputDiversity <- list("1" = numeric(0), "2" = numeric(0), "3" = numeric(0),
                        "4" = numeric(0), "5" = numeric(0), "6" = numeric(0),
                        "7" = numeric(0), "8" = numeric(0), "9" = numeric(0), 
                        "10" = numeric(0))
for (i in 1:n) {
  listOutput <- list()
  listOutput$"50" <- TS_Algorithm(1, 50, nodes, countIDs, method, randomize)
  listOutput$"100" <- TS_Algorithm(1, 100, nodes, countIDs, method, randomize)
  listOutput$"150" <- TS_Algorithm(1, 150, nodes, countIDs, method, randomize)
  listOutput$"200" <- TS_Algorithm(1, 200, nodes, countIDs, method, randomize)
  listOutput$"250" <- TS_Algorithm(1, 250, nodes, countIDs, method, randomize)
  listOutput$"300" <- TS_Algorithm(1, 300, nodes, countIDs, method, randomize)
  listOutput$"350" <- TS_Algorithm(1, 350, nodes, countIDs, method, randomize)
  listOutput$"400" <- TS_Algorithm(1, 400, nodes, countIDs, method, randomize)
#}  
  evalOutput <- list()
  evalOutput$"50" <- Evaluate_TS(listOutput$"50", nodes, countIDs)
  evalOutput$"100" <- Evaluate_TS(listOutput$"100", nodes, countIDs)
  evalOutput$"150" <- Evaluate_TS(listOutput$"150", nodes, countIDs)
  evalOutput$"200" <- Evaluate_TS(listOutput$"200", nodes, countIDs)
  evalOutput$"250" <- Evaluate_TS(listOutput$"250", nodes, countIDs)
  evalOutput$"300" <- Evaluate_TS(listOutput$"300", nodes, countIDs)
  evalOutput$"350" <- Evaluate_TS(listOutput$"350", nodes, countIDs)
  evalOutput$"400" <- Evaluate_TS(listOutput$"400", nodes, countIDs)
  
  
  for (level in 3:10) {
    diversity <- numeric(0)
    for (number in names(evalOutput)) {
      diversity <- c(diversity, length(evalOutput[[number]][[level]]))
    }
    outputDiversity[[level]] <- rbind(outputDiversity[[level]], diversity)
  }
}

method <- "balance"
outputBalance <- list("1" = numeric(0), "2" = numeric(0), "3" = numeric(0),
                      "4" = numeric(0), "5" = numeric(0), "6" = numeric(0),
                      "7" = numeric(0), "8" = numeric(0), "9" = numeric(0), 
                      "10" = numeric(0))
for (i in 1:n) {
  listOutput <- list()
  listOutput$"50" <- TS_Algorithm(1, 50, nodes, countIDs, method, randomize)
  listOutput$"100" <- TS_Algorithm(1, 100, nodes, countIDs, method, randomize)
  listOutput$"150" <- TS_Algorithm(1, 150, nodes, countIDs, method, randomize)
  listOutput$"200" <- TS_Algorithm(1, 200, nodes, countIDs, method, randomize)
  listOutput$"250" <- TS_Algorithm(1, 250, nodes, countIDs, method, randomize)
  listOutput$"300" <- TS_Algorithm(1, 300, nodes, countIDs, method, randomize)
  listOutput$"350" <- TS_Algorithm(1, 350, nodes, countIDs, method, randomize)
  listOutput$"400" <- TS_Algorithm(1, 400, nodes, countIDs, method, randomize)
  
  evalOutput <- list()
  evalOutput$"50" <- Evaluate_TS(listOutput$"50", nodes, countIDs)
  evalOutput$"100" <- Evaluate_TS(listOutput$"100", nodes, countIDs)
  evalOutput$"150" <- Evaluate_TS(listOutput$"150", nodes, countIDs)
  evalOutput$"200" <- Evaluate_TS(listOutput$"200", nodes, countIDs)
  evalOutput$"250" <- Evaluate_TS(listOutput$"250", nodes, countIDs)
  evalOutput$"300" <- Evaluate_TS(listOutput$"300", nodes, countIDs)
  evalOutput$"350" <- Evaluate_TS(listOutput$"350", nodes, countIDs)
  evalOutput$"400" <- Evaluate_TS(listOutput$"400", nodes, countIDs)
  
  
  for (level in 3:10) {
    diversity <- numeric(0)
    for (number in names(evalOutput)) {
      diversity <- c(diversity, length(evalOutput[[number]][[level]]))
    }
    outputBalance[[level]] <- rbind(outputBalance[[level]], diversity)
  }
}


outputRandom <- list("1" = numeric(0), "2" = numeric(0), "3" = numeric(0),
                     "4" = numeric(0), "5" = numeric(0), "6" = numeric(0),
                     "7" = numeric(0), "8" = numeric(0), "9" = numeric(0), 
                     "10" = numeric(0))
for (i in 1:n) {
  listRandom <- list()
  listRandom$"50" <- RandomSampling(idsFile, nodes, 50)
  listRandom$"100" <- RandomSampling(idsFile, nodes, 100)
  listRandom$"150" <- RandomSampling(idsFile, nodes, 150)
  listRandom$"200" <- RandomSampling(idsFile, nodes, 200)
  listRandom$"250" <- RandomSampling(idsFile, nodes, 250)
  listRandom$"300" <- RandomSampling(idsFile, nodes, 300)
  listRandom$"350" <- RandomSampling(idsFile, nodes, 350)
  listRandom$"400" <- RandomSampling(idsFile, nodes, 400)
  
  evalRandom <- list()
  evalRandom$"50" <- Evaluate_TS(listRandom$"50", nodes, countIDs)
  evalRandom$"100" <- Evaluate_TS(listRandom$"100", nodes, countIDs)
  evalRandom$"150" <- Evaluate_TS(listRandom$"150", nodes, countIDs)
  evalRandom$"200" <- Evaluate_TS(listRandom$"200", nodes, countIDs)
  evalRandom$"250" <- Evaluate_TS(listRandom$"250", nodes, countIDs)
  evalRandom$"300" <- Evaluate_TS(listRandom$"300", nodes, countIDs)
  evalRandom$"350" <- Evaluate_TS(listRandom$"350", nodes, countIDs)
  evalRandom$"400" <- Evaluate_TS(listRandom$"400", nodes, countIDs)

  for (level in 3:10) {
    diversity <- numeric(0)
    for (number in names(evalRandom)) {
      diversity <- c(diversity, length(evalRandom[[number]][[level]]))
    }
    outputRandom[[level]] <- rbind(outputRandom[[level]], diversity)
  }
}


totalTaxa <- numeric(0)
children <- 1
for (i in 1:10) {
  taxon <- as.integer(children)
  children <- nodes$id[is.element(nodes$parent, taxon) & 
                       !is.element(nodes$id, taxon)]
  children <- intersect(children, names(countIDs))
  totalTaxa <- c(totalTaxa, length(children))
}


save.image("Eval_final.RData")

confidence <- .995  # 99% = .995, 95% = .975

diversityMeans <- list()
diversityCI <- list()

balanceMeans <- list()
balanceCI <- list()

randomMeans <- list()
randomCI <- list()

for (level in 3:10) {
  diversityMeans[[level]] <- colMeans(outputDiversity[[level]])
  balanceMeans[[level]] <- colMeans(outputBalance[[level]])
  randomMeans[[level]] <- colMeans(outputRandom[[level]])

  diversityCI[[level]] <- apply(outputDiversity[[level]], 2, sd)
  balanceCI[[level]] <- apply(outputBalance[[level]], 2, sd)
  randomCI[[level]] <- apply(outputRandom[[level]], 2, sd)

  diversityCI[[level]] <- diversityCI[[level]]/sqrt(n)
  balanceCI[[level]] <- balanceCI[[level]]/sqrt(n)
  randomCI[[level]] <- randomCI[[level]]/sqrt(n)

  diversityCI[[level]] <- qt(confidence, df = n-1) * diversityCI[[level]]
  balanceCI[[level]] <- qt(confidence, df = n-1) * balanceCI[[level]]
  randomCI[[level]] <- qt(confidence, df = n-1) * randomCI[[level]]
}


for (level in 3:10) {
  df <- data.frame(x, diversityMeans = diversityMeans[[level]], 
                      balanceMeans = balanceMeans[[level]],
                      randomMeans = randomMeans[[level]],
                      totalTaxa = rep(totalTaxa[level], 8))
  
  imageName <- paste0("level", level, ".pdf")
  pdf(imageName)
  print(ggplot(df, aes(x)) + 
          geom_point(aes(y=diversityMeans, colour="TSdiversity")) +
          geom_line(aes(y=diversityMeans, colour="TSdiversity")) + 
          geom_errorbar(aes(ymin=diversityMeans - diversityCI[[level]],
                            ymax=diversityMeans + diversityCI[[level]],
                            colour = "TSdiversity"), width=1) +
          geom_point(aes(y=balanceMeans, colour="TSbalance")) +
          geom_line(aes(y=balanceMeans, colour="TSbalance")) + 
          geom_errorbar(aes(ymin=balanceMeans - balanceCI[[level]],
                            ymax=balanceMeans + balanceCI[[level]],
                            colour = "TSbalance"), width=1) +
          geom_point(aes(y=randomMeans, colour="RS")) +
          geom_line(aes(y=randomMeans, colour="RS")) +
          geom_errorbar(aes(ymin=randomMeans - randomCI[[level]],
                            ymax=randomMeans + randomCI[[level]],
                            colour = "RS"), width=1) +
          geom_point(aes(y=totalTaxa, colour="max")) +
          geom_line(aes(y=totalTaxa, colour="max")) +
          xlab("m") + ylab(paste0("# taxa (level = ", level, ")")) +
          scale_color_manual("Method",
                             values = c("TSdiversity" = "orange",
                                        "TSbalance" = "red",
                                        "RS" = "blue",
                                        "max" = "black")))
  dev.off()
}
