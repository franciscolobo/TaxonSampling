
#' Retrieve taxonomy and phylogeny datasets
#'
#' Downloads and extracts datasets from NCBI Taxonomy,
#' Open Tree Taxonomy (OTT), and the OpenTree synthetic phylogeny.
#'
#' All three datasets are downloaded by default.
#'
#' Compressed archives are extracted automatically and deleted
#' after successful extraction. Uncompressed files are preserved.
#'
#' No taxonomic or phylogenetic processing is performed.
#'
#' @param target.dir Directory where datasets will be saved.
#' @param datasets Character vector containing any combination of
#'   "ncbi", "ott", and "opentree". Defaults to all three.
#' @param method Download method passed to utils::download.file().
#' @param unzip Unzip method passed to utils::unzip().
#' @param overwrite Logical. Replace existing datasets if TRUE.
#' @param ncbi.url URL of the NCBI Taxonomy archive.
#' @param ott.url URL of the Open Tree Taxonomy archive.
#' @param opentree.url URL of the OpenTree synthetic phylogeny.
#' @param quiet Logical. Suppress progress messages.
#' @param ... Additional arguments passed to download.file().
#'
#' @return Invisibly returns a data.frame with dataset names,
#'   output directories, source URLs, and download status.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' retrieve_data_files("data_files")
#' retrieve_data_files("data_files", datasets = "ncbi")
#' retrieve_data_files("data_files", overwrite = TRUE)
#' }

retrieve_taxonomy_files <- function(
    target.dir,
    datasets = c("ncbi", "ott", "opentree"),
    method = "auto",
    unzip = getOption("unzip"),
    overwrite = FALSE,
    ncbi.url = paste0(
      "https://ftp.ncbi.nlm.nih.gov/",
      "pub/taxonomy/taxdump.tar.gz"
    ),
    ott.url = paste0(
      "https://files.opentreeoflife.org/",
      "ott/ott3.7/ott3.7.tgz"
    ),
    opentree.url = paste0(
      "https://files.opentreeoflife.org/",
      "synthesis/current/output/",
      "labelled_supertree/",
      "labelled_supertree_ottnames.tre"
    ),
    quiet = FALSE,
    ...
) {
  
  # ============================================================
  # 1. Validate arguments
  # ============================================================
  
  if (!is.character(target.dir) ||
      length(target.dir) != 1L ||
      is.na(target.dir) ||
      !nzchar(target.dir)) {
    stop("target.dir must be a non-empty character string.")
  }
  
  allowed <- c("ncbi", "ott", "opentree")
  
  if (!is.character(datasets) ||
      !length(datasets) ||
      anyNA(datasets) ||
      any(!datasets %in% allowed)) {
    stop(
      "datasets must contain one or more of: ",
      paste(allowed, collapse = ", ")
    )
  }
  
  datasets <- unique(datasets)
  
  for (x in list(overwrite = overwrite, quiet = quiet)) {
    if (!is.logical(x) || length(x) != 1L || is.na(x)) {
      stop("overwrite and quiet must be TRUE or FALSE.")
    }
  }
  
  if (!dir.exists(target.dir)) {
    dir.create(
      target.dir,
      recursive = TRUE,
      showWarnings = FALSE
    )
  }
  
  if (!dir.exists(target.dir)) {
    stop("Could not create target.dir: ", target.dir)
  }
  
  target.dir <- normalizePath(
    target.dir,
    mustWork = TRUE
  )
  
  # ============================================================
  # 2. Dataset sources
  # ============================================================
  
  sources <- list(
    ncbi = ncbi.url,
    ott = ott.url,
    opentree = opentree.url
  )
  
  # ============================================================
  # 3. Helper: identify archive format
  # ============================================================
  
  archive_type <- function(filename) {
    
    filename <- tolower(filename)
    
    if (grepl("\\.zip$", filename)) {
      return("zip")
    }
    
    if (grepl("\\.(tar\\.gz|tgz)$", filename)) {
      return("tar.gz")
    }
    
    if (grepl("\\.tar$", filename)) {
      return("tar")
    }
    
    return("none")
  }
  
  # ============================================================
  # 4. Helper: validate archive member paths
  # ============================================================
  
  validate_members <- function(members) {
    
    members <- gsub("\\\\", "/", members)
    
    unsafe <- grepl(
      "^/|^[A-Za-z]:|(^|/)\\.\\.(/|$)",
      members
    )
    
    if (any(unsafe)) {
      stop(
        "Archive contains unsafe paths: ",
        paste(head(members[unsafe], 5L), collapse = ", ")
      )
    }
    
    invisible(TRUE)
  }
  
  # ============================================================
  # 5. Helper: extract compressed files
  # ============================================================
  
  extract_archive <- function(
    archive,
    destination,
    type
  ) {
    
    if (type == "zip") {
      
      members <- utils::unzip(
        archive,
        list = TRUE
      )$Name
      
      validate_members(members)
      
      utils::unzip(
        archive,
        exdir = destination,
        unzip = unzip
      )
      
    } else if (type %in% c("tar.gz", "tar")) {
      
      members <- utils::untar(
        archive,
        list = TRUE
      )
      
      validate_members(members)
      
      utils::untar(
        archive,
        exdir = destination
      )
      
    } else {
      
      stop("Unsupported archive type: ", type)
    }
    
    # Remove macOS metadata directories if present.
    macos <- file.path(destination, "__MACOSX")
    
    if (dir.exists(macos)) {
      unlink(
        macos,
        recursive = TRUE,
        force = TRUE
      )
    }
    
    # Verify that extraction produced files.
    extracted <- list.files(
      destination,
      recursive = TRUE,
      full.names = TRUE,
      all.files = TRUE,
      no.. = TRUE
    )
    
    extracted <- extracted[
      !dir.exists(extracted)
    ]
    
    if (!length(extracted)) {
      stop(
        "Archive extraction produced no files: ",
        archive
      )
    }
    
    invisible(extracted)
  }
  
  # ============================================================
  # 6. Helper: download and install one dataset
  # ============================================================
  
  download_dataset <- function(dataset) {
    
    url <- sources[[dataset]]
    
    if (!is.character(url) ||
        length(url) != 1L ||
        is.na(url) ||
        !grepl("^https?://", url)) {
      stop("Invalid URL for dataset: ", dataset)
    }
    
    filename <- basename(
      sub("[?#].*$", "", url)
    )
    
    if (!nzchar(filename) ||
        filename %in% c(".", "..")) {
      stop("Invalid filename for dataset: ", dataset)
    }
    
    destination <- file.path(
      target.dir,
      dataset
    )
    
    # ----------------------------------------------------------
    # Skip existing datasets
    # ----------------------------------------------------------
    
    if (dir.exists(destination) && !overwrite) {
      
      if (!quiet) {
        message(
          "[", dataset, "] Already exists. Skipping."
        )
      }
      
      return(data.frame(
        dataset = dataset,
        directory = destination,
        url = url,
        status = "skipped",
        stringsAsFactors = FALSE
      ))
    }
    
    # ----------------------------------------------------------
    # Create staging directory
    # ----------------------------------------------------------
    
    staging <- tempfile(
      pattern = paste0(dataset, "_staging_"),
      tmpdir = target.dir
    )
    
    dir.create(staging)
    
    on.exit(
      unlink(
        staging,
        recursive = TRUE,
        force = TRUE
      ),
      add = TRUE
    )
    
    archive <- file.path(
      staging,
      filename
    )
    
    output.dir <- file.path(
      staging,
      "extracted"
    )
    
    dir.create(output.dir)
    
    # ----------------------------------------------------------
    # Download original file
    # ----------------------------------------------------------
    
    if (!quiet) {
      message(
        "[", dataset, "] Downloading: ", url
      )
    }
    
    status <- tryCatch(
      utils::download.file(
        url = url,
        destfile = archive,
        method = method,
        mode = "wb",
        quiet = quiet,
        cacheOK = FALSE,
        ...
      ),
      error = function(e) {
        stop(
          "Download failed for ", dataset, ": ",
          conditionMessage(e)
        )
      }
    )
    
    if (!identical(as.integer(status), 0L) ||
        !file.exists(archive) ||
        is.na(file.size(archive)) ||
        file.size(archive) == 0) {
      stop("Invalid or incomplete download: ", dataset)
    }
    
    # ----------------------------------------------------------
    # Extract archives or retain uncompressed files
    # ----------------------------------------------------------
    
    type <- archive_type(filename)
    
    if (type == "none") {
      
      # Move the original uncompressed file unchanged.
      if (!file.rename(
        archive,
        file.path(output.dir, filename)
      )) {
        stop("Could not move downloaded file: ", dataset)
      }
      
    } else {
      
      if (!quiet) {
        message(
          "[", dataset, "] Extracting archive..."
        )
      }
      
      extract_archive(
        archive = archive,
        destination = output.dir,
        type = type
      )
      
      # Archive is no longer needed.
      unlink(archive)
      
      if (!quiet) {
        message(
          "[", dataset, "] Extraction complete."
        )
      }
    }
    
    # ----------------------------------------------------------
    # Install completed dataset
    # ----------------------------------------------------------
    
    backup <- NULL
    
    if (dir.exists(destination)) {
      
      backup <- tempfile(
        pattern = paste0(dataset, "_backup_"),
        tmpdir = target.dir
      )
      
      if (!file.rename(destination, backup)) {
        stop(
          "Could not back up existing dataset: ",
          dataset
        )
      }
    }
    
    installed <- file.rename(
      output.dir,
      destination
    )
    
    if (!installed) {
      
      # Restore previous dataset.
      if (!is.null(backup)) {
        if (!file.rename(backup, destination)) {
          warning(
            "Could not restore backup at: ",
            backup
          )
        }
      }
      
      stop(
        "Could not install dataset: ",
        dataset
      )
    }
    
    # Delete previous version after successful installation.
    if (!is.null(backup)) {
      unlink(
        backup,
        recursive = TRUE,
        force = TRUE
      )
    }
    
    if (!quiet) {
      message(
        "[", dataset, "] Saved to: ",
        destination
      )
    }
    
    data.frame(
      dataset = dataset,
      directory = destination,
      url = url,
      status = "downloaded",
      stringsAsFactors = FALSE
    )
  }
  
  # ============================================================
  # 7. Download requested datasets
  # ============================================================
  
  results <- lapply(
    datasets,
    download_dataset
  )
  
  results <- do.call(
    rbind,
    results
  )
  
  rownames(results) <- NULL
  
  # ============================================================
  # 8. Summary
  # ============================================================
  
  if (!quiet) {
    
    message("\nDataset retrieval completed.")
    
    message(
      "Downloaded: ",
      sum(results$status == "downloaded")
    )
    
    message(
      "Skipped: ",
      sum(results$status == "skipped")
    )
    
    message(
      "Directory: ",
      target.dir
    )
  }
  
  invisible(results)
}
