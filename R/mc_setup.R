mc_setup <- function(ncpus){
  if(ncpus > parallel::detectCores() - 1){
    ncpus <- parallel::detectCores() - 1
    warning("Note: Only ", parallel::detectCores(), " nodes available. Using ",
            ncpus, " for counting species.")
  }

  if (ncpus > 1){
    if(.Platform$OS.type == "windows"){
      cl <- parallel::makeCluster(ncpus, setup_strategy = "sequential")
      parallel::clusterExport(cl, "nodes")
    } else {
      cl <- ncpus
    }
  } else {
    cl <- 1
  }
  return(cl)
}
