#evaluate_TS <- function(x) {
  # Collect metrics about the taxon diversity of a sample (outputIDs).
  #
  # Args:
  #   taxlist: _taxonsampling_ list object returned by [run_TS()]
  # Returns:
  #   listTaxon: (char) children taxa selected at least once per level.
  ow <- options("warn")
  options(warn = -1)
  selectedIDs <-  get_taxonomy_counts(ids_df = data.frame(taxID = x$outputIDs,
                                                          seqID = NA),
                                      nodes = x$nodes, verbose = FALSE)
  options(ow)
  listTaxon   <- list()
  listBias    <- list()

  # Find the sub-taxa (children nodes) of the current taxon
  children <- 1
  while (length(children) > 0) {
    taxon    <- as.integer(children)
    children <- selectedIDs$nodes$id[(selectedIDs$nodes$parent %in% taxon) & !(selectedIDs$nodes$id %in% taxon)]
    children <- intersect(children, names(selectedIDs$countIDs))

    selectedChildren <- intersect(children, names(selectedIDs$countIDs))
    listTaxon        <- c(listTaxon, list(selectedChildren))
  }

  return(listTaxon)
}


library(dplyr)

#Number of reps
n <- 10

# preprocess stuff:
taxlist <- get_counts(taxonomy_path = "data_files/taxdump/",
                      ids_file      = "data_files/metadata/TaxID2SeqID.txt")

# Test parameters
rand <- c("yes", "no", "after_first_round")
meth <- c("diversity", "balanced")
samp <- c("agnostic", "known_species")
repl <- c(TRUE, FALSE)
m    <- 50 * (1:8)
pars <- expand.grid(rand=rand, meth=meth, samp=samp, repl=repl, m = m,
                    stringsAsFactors = FALSE)
pars <- pars[-which(pars$meth == "balanced" & pars$rand == "after_first_round"), ]

output <- vector("list", nrow(pars))
for (k in 1:nrow(pars)){
  cat(sprintf("\nTrying %02d/%02d: [%s]", k, nrow(pars), paste(pars[k, ], collapse = ",")))
  output[[k]] <- list(pars = pars[k, ])

  out <- vector("list", n)
  names(out) <- paste0("Rep", 1:n)
  for (i in 1:n){
    out[[i]]$tl <- run_TS(taxlist          = taxlist,
                          taxon            = 40674,
                          m                = pars$m[k],
                          seq_file         = NULL,
                          out_file         = NULL,
                          method           = pars$meth[k],
                          randomize        = pars$rand[k],
                          replacement      = pars$repl[k],
                          ignoreIDs        = 10090,
                          requireIDs       = 9606,
                          sampling         = pars$samp[k],
                          verbose          = FALSE)

    #out[[i]]$perf <- evaluate_TS(out[[i]]$tl)
    out[[i]]$tl$nodes  <- NULL
    out[[i]]$tl$spp_df <- NULL
  }

  output[[k]]$out <- out
}

saveRDS(output, "./results/test.rds")




