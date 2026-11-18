#' Download and build the available TaxonSampling reference databases
#' @param path Parent directory for source files and databases.
#' @param sources Sources to retrieve: ncbi, ott, synthetic (or opentree).
#' @param download Download sources before building (FALSE uses existing files).
#' @param build_cross_identity Build identifier mapping when possible.
#' @param overwrite Replace existing sources/databases.
#' @param include_suppressed Include OTT suppressed taxa in species counts.
#' @param verbose Display progress messages.
#' @param progress_callback Optional function receiving a human-readable status message.
#' @return Invisibly returns a data.frame describing each database build.
#' @export
setup_taxonsampling <- function(path = "taxonomy_data",
                                sources = c("ncbi", "ott", "synthetic"),
                                download = TRUE,
                                build_cross_identity = TRUE,
                                overwrite = FALSE,
                                include_suppressed = TRUE,
                                verbose = TRUE,
                                progress_callback = NULL) {
  notify <- function(text) {
    if (is.function(progress_callback)) progress_callback(text)
  }
  allowed <- c("ncbi", "ott", "synthetic")
  sources[sources == "opentree"] <- "synthetic"
  if (!length(sources) || anyNA(sources) || any(!sources %in% allowed))
    stop("sources must be ncbi, ott, or synthetic")
  sources <- unique(sources)
  raw_dir <- file.path(path, "raw")
  db_dir <- file.path(path, "taxonomy_db")
  dir.create(raw_dir, recursive=TRUE, showWarnings=FALSE)
  dir.create(db_dir, recursive=TRUE, showWarnings=FALSE)
  if (download) {
    notify("Starting taxonomy downloads")
    retrieve_taxonomy_files(raw_dir,
      datasets=ifelse(sources=="synthetic", "opentree", sources),
      overwrite=overwrite, quiet=!verbose, continue_on_error=TRUE,
      progress_callback=progress_callback)
  }
  input_dirs <- c(ncbi=file.path(raw_dir,"ncbi"),
                  ott=file.path(raw_dir,"ott"),
                  synthetic=file.path(raw_dir,"opentree"))
  # A synthetic tree may be nested in the downloaded directory.
  synthetic_path <- function(d) {
    if (!dir.exists(d)) return(NA_character_)
    matches <- list.files(d, pattern="\\.(tre|tree|nwk|newick)$",
                          recursive=TRUE, full.names=TRUE, ignore.case=TRUE)
    if (length(matches)!=1L) stop("Expected exactly one synthetic Newick file in ", d)
    matches[[1L]]
  }
  result <- list()
  for (src in sources) {
    db <- file.path(db_dir, paste0(src,".duckdb"))
    if (file.exists(db) && !overwrite) {
      notify(paste0("Using existing database: ", src))
      result[[src]] <- data.frame(database=src, status="existing", path=db)
      next
    }
    input <- input_dirs[[src]]
    if (!dir.exists(input)) {
      warning("Skipping ", src, ": downloaded files unavailable", call.=FALSE)
      result[[src]] <- data.frame(database=src, status="missing_source", path=db)
      next
    }
    if (src=="synthetic") {
      input <- tryCatch(synthetic_path(input), error=function(e) {
        warning(conditionMessage(e), call.=FALSE); NA_character_
      })
      if (is.na(input)) {
        result[[src]] <- data.frame(database=src, status="missing_tree", path=db)
        next
      }
    }
    notify(paste0("Creating database: ", src, ".duckdb"))
    status <- tryCatch({
      build_taxonomy_db(taxonomy_path=input, taxonomy=src, output_dir=db_dir,
                        ott_path=if (src=="synthetic" && dir.exists(input_dirs[["ott"]]))
                          input_dirs[["ott"]] else NULL,
                        include_suppressed=include_suppressed,
                        overwrite=overwrite, verbose=verbose)
      "built"
    }, error=function(e) {
      warning("Failed to build ", src, ": ", conditionMessage(e), call.=FALSE)
      "failed"
    })
    notify(paste0("Database ", src, ": ", status))
    result[[src]] <- data.frame(database=src, status=status, path=db)
  }
  cross_db <- file.path(db_dir, "cross_identity.duckdb")
  if (build_cross_identity) {
    if (file.exists(cross_db) && !overwrite) {
      cross_status <- "existing"
      notify("Using existing database: cross_identity")
      if (any(vapply(result, function(x) x$status=="built", logical(1)))) {
        warning("Cross-identity database was not rebuilt despite new taxonomy builds; ",
                "consider overwrite=TRUE to synchronize versions.", call.=FALSE)
      }
    } else {
      notify("Creating database: cross_identity.duckdb")
      cross_status <- tryCatch({
        build_cross_identity_db(db_dir, output_path=cross_db, overwrite=overwrite,
                                verbose=verbose)
        "built"
      }, error=function(e) {
        warning("Cross-identity database unavailable: ", conditionMessage(e), call.=FALSE)
        "unavailable"
      })
    }
    notify(paste0("Database cross_identity: ", cross_status))
    result[["cross_identity"]] <- data.frame(database="cross_identity",
                                               status=cross_status, path=cross_db)
  }
  out <- do.call(rbind, result)
  rownames(out) <- NULL
  if (verbose) print(out, row.names=FALSE)
  notify("Taxonomy database setup finished")
  invisible(out)
}

#' Inspect TaxonSampling reference database availability
#' @param path Parent directory used by setup_taxonsampling().
#' @return Data frame with availability and selected metadata.
#' @export
taxonomy_status <- function(path = "taxonomy_data") {
  db_dir <- file.path(path, "taxonomy_db")
  names <- c("ncbi", "ott", "synthetic", "cross_identity")
  result <- lapply(names, function(nm) {
    file <- file.path(db_dir, paste0(nm, ".duckdb"))
    if (!file.exists(file)) return(data.frame(database=nm, available=FALSE,
                                              build_date=NA_character_, path=file))
    build_date <- tryCatch({
      con <- DBI::dbConnect(duckdb::duckdb(), dbdir=file, read_only=TRUE)
      on.exit(DBI::dbDisconnect(con, shutdown=TRUE), add=TRUE)
      x <- DBI::dbGetQuery(con, "SELECT value FROM metadata WHERE key='build_date'")
      if (nrow(x)) as.character(x$value[[1L]]) else NA_character_
    }, error=function(e) NA_character_)
    data.frame(database=nm, available=TRUE, build_date=build_date, path=file)
  })
  do.call(rbind, result)
}
