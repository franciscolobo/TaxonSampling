#' Build a reusable TaxonSampling scaffold database
#'
#' Reads a complete taxonomy or Newick tree, computes subtree species and tip
#' counts independently of any input sequence identifiers, and writes a DuckDB
#' database. NCBI/OTT cross-references are stored when available, but NEVER
#' determine which nodes contribute to background counts.
#'
#' @param taxonomy_path Directory containing extracted taxonomy files, or a
#'   Newick file for `custom` (also accepted for `synthetic`).
#' @param taxonomy One of "ncbi", "ott", "synthetic", "custom".
#' @param output_dir Directory in which to create the database.
#' @param db_name Optional output database filename. Defaults to taxonomy.duckdb.
#' @param ott_path Optional OTT taxonomy directory/file, needed to annotate
#'   synthetic-tree species ranks. Without it, synthetic species_count is NA.
#' @param include_suppressed Whether OTT flagged/suppressed taxa count toward
#'   species_count. This only affects the OTT taxonomy (and synthetic when OTT
#'   metadata are provided). Default TRUE retains the complete background.
#' @param overwrite Replace an existing output database.
#' @param verbose Print progress messages.
#' @return Invisibly returns the database path.
#' @export
build_taxonomy_db <- function(taxonomy_path,
                              taxonomy = c("ncbi", "ott", "synthetic", "custom"),
                              output_dir = "taxonomy_db",
                              db_name = NULL,
                              ott_path = NULL,
                              include_suppressed = TRUE,
                              overwrite = FALSE,
                              verbose = TRUE) {
  taxonomy <- match.arg(taxonomy)
  for (pkg in c("data.table", "DBI", "duckdb")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package '", pkg, "' is required.", call. = FALSE)
    }
  }
  if (taxonomy %in% c("synthetic", "custom") &&
      !requireNamespace("ape", quietly = TRUE)) {
    stop("Package 'ape' is required for Newick trees.", call. = FALSE)
  }
  if (!is.character(taxonomy_path) || length(taxonomy_path) != 1L ||
      !file.exists(taxonomy_path)) stop("taxonomy_path does not exist.")
  if (!is.logical(include_suppressed) || length(include_suppressed) != 1L ||
      is.na(include_suppressed)) stop("include_suppressed must be TRUE or FALSE.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (is.null(db_name)) db_name <- paste0(taxonomy, ".duckdb")
  if (basename(db_name) != db_name || !grepl("\\.duckdb$", db_name)) {
    stop("db_name must be a filename ending in .duckdb")
  }
  destination <- normalizePath(file.path(output_dir, db_name), mustWork = FALSE)
  if (file.exists(destination) && !overwrite) {
    stop("Database already exists: ", destination, "; set overwrite=TRUE")
  }
  say <- function(...) if (isTRUE(verbose)) message(...)
  dt <- data.table::as.data.table
  trim <- function(x) trimws(as.character(x))
  
  # Discover one source file in the supplied extracted directory, even when
  # the archive includes an enclosing versioned directory.
  locate <- function(path, filename) {
    if (file.exists(path) && !dir.exists(path)) {
      if (basename(path) == filename) return(normalizePath(path))
      stop("Expected file named '", filename, "', got: ", path)
    }
    matches <- list.files(path, pattern = paste0("^", gsub("\\.", "\\\\.", filename), "$"),
                          recursive = TRUE, full.names = TRUE)
    if (length(matches) != 1L) {
      stop("Expected exactly one '", filename, "' under ", path,
           "; found ", length(matches), ".")
    }
    normalizePath(matches)
  }
  read_ott <- function(path) {
    f <- locate(path, "taxonomy.tsv")
    say("Reading OTT taxonomy: ", f)
    # OTT uses tab-pipe-tab separators, with a trailing pipe on each row.
    x <- data.table::fread(f, sep = "|", quote = "", fill = TRUE,
                           showProgress = verbose, colClasses = "character")
    data.table::setnames(x, trim(names(x)))
    wanted <- c("uid", "parent_uid", "name", "rank")
    if (!all(wanted %in% names(x))) {
      stop("OTT taxonomy.tsv must contain: ", paste(wanted, collapse = ", "),
           ". Found: ", paste(names(x), collapse = ", "))
    }
    out <- data.table::data.table(
      id = trim(x[["uid"]]), parent = trim(x[["parent_uid"]]),
      name = trim(x[["name"]]), rank = tolower(trim(x[["rank"]])),
      flags = if ("flags" %in% names(x)) trim(x[["flags"]]) else NA_character_,
      sourceinfo = if ("sourceinfo" %in% names(x)) trim(x[["sourceinfo"]]) else NA_character_
    )
    out[parent %in% c("", "0", "NA"), parent := NA_character_]
    out
  }
  read_ncbi <- function(path) {
    nf <- locate(path, "nodes.dmp")
    sf <- locate(path, "names.dmp")
    say("Reading NCBI nodes: ", nf)
    n <- data.table::fread(nf, sep = "|", quote = "", header = FALSE,
                           select = 1:3, colClasses = "character", showProgress = verbose)
    data.table::setnames(n, c("id", "parent", "rank"))
    for (j in names(n)) data.table::set(n, j = j, value = trim(n[[j]]))
    say("Reading NCBI scientific names: ", sf)
    s <- data.table::fread(sf, sep = "|", quote = "", header = FALSE,
                           select = c(1L, 2L, 4L), colClasses = "character",
                           showProgress = verbose)
    data.table::setnames(s, c("id", "name", "class"))
    for (j in names(s)) data.table::set(s, j = j, value = trim(s[[j]]))
    s <- s[class == "scientific name", .(id, name)]
    if (anyDuplicated(s$id)) stop("Duplicate NCBI scientific names for a TaxID.")
    out <- merge(n, s, by = "id", all.x = TRUE, sort = FALSE)
    out[, `:=`(flags = NA_character_, sourceinfo = NA_character_)]
    out
  }
  # OpenTree labels may be unquoted (Name_ott123) or quoted in Newick
  # ('Name with spaces ott123'). Some ape versions retain a closing quote
  # after parsing, so accept it without changing the original node label.
  extract_ott_id <- function(labels) {
    labels <- trimws(as.character(labels))
    # ape may preserve a trailing quote while removing the leading one.
    # Normalize only for identifier extraction; keep original names intact.
    labels <- sub("^'+", "", labels)
    labels <- sub("'+$", "", labels)
    result <- rep(NA_character_, length(labels))
    valid <- which(!is.na(labels))
    if (!length(valid)) return(result)
    # Require a separator before ott, so synthetic MRCA labels never match.
    matches <- regexec("(?:_|[[:space:]])ott([0-9]+)$", labels[valid], perl = TRUE)
    pieces <- regmatches(labels[valid], matches)
    has_id <- lengths(pieces) >= 2L
    if (any(has_id)) {
      result[valid[has_id]] <- vapply(pieces[has_id], `[[`, character(1), 2L)
    }
    result
  }
  newick_warnings <- character()
  read_newick <- function(path) {
    if (dir.exists(path)) {
      candidates <- list.files(path, recursive = TRUE, full.names = TRUE,
                               pattern = "\\.(tre|tree|nwk|newick)$", ignore.case = TRUE)
      if (length(candidates) != 1L) {
        stop("Expected one Newick file in directory; found ", length(candidates),
             ". Supply the file path explicitly.")
      }
      path <- candidates
    }
    say("Reading Newick tree: ", path)
    # ape may warn about numeric coercion for some OpenTree labels. Preserve
    # warnings in metadata; do not suppress unrelated warnings.
    tree <- withCallingHandlers(
      ape::read.tree(path),
      warning = function(w) {
        if (identical(conditionMessage(w), "NAs introduced by coercion") &&
            grepl("FOO\\(lab\\[k\\]\\)",
                  paste(deparse(conditionCall(w)), collapse = ""))) {
          newick_warnings <<- c(newick_warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      }
    )
    if (inherits(tree, "multiPhylo")) stop("Expected a single Newick tree.")
    nt <- length(tree$tip.label)
    ni <- tree$Nnode
    labels <- c(tree$tip.label, rep(NA_character_, ni))
    if (!is.null(tree$node.label)) {
      if (length(tree$node.label) != ni) stop("Unexpected internal label count.")
      labels[nt + seq_len(ni)] <- tree$node.label
    }
    # The OTT suffix is a stable taxon ID, even when the taxon name is repeated.
    ott <- extract_ott_id(labels)
    # mrcaott... labels are synthetic topology nodes, not OTT taxon IDs.
    ids <- ifelse(!is.na(ott), paste0("ott:", ott), NA_character_)
    labelled <- !is.na(labels) & nzchar(labels)
    ids[is.na(ids) & labelled] <- paste0("label:", labels[is.na(ids) & labelled])
    ids[is.na(ids)] <- paste0("node:", which(is.na(ids)))
    if (anyDuplicated(ids)) {
      stop("Newick labels are not unique. Supply unique node labels or OTT IDs.")
    }
    parent_idx <- rep(NA_integer_, length(labels))
    parent_idx[tree$edge[, 2L]] <- tree$edge[, 1L]
    node_type <- ifelse(seq_along(labels) <= nt, "tip", "internal")
    node_type[grepl("^mrcaott[0-9]+ott[0-9]+$", labels) & !is.na(labels)] <- "synthetic_mrca"
    node_type[!labelled & seq_along(labels) > nt] <- "unlabeled"
    label_status <- ifelse(!labelled, "missing",
                           ifelse(!is.na(ott), "ott_taxon",
                                  ifelse(node_type == "synthetic_mrca",
                                         "synthetic_mrca", "unrecognized")))
    data.table::data.table(
      id = ids, parent = ifelse(is.na(parent_idx), NA_character_, ids[parent_idx]),
      name = labels, rank = NA_character_, flags = NA_character_,
      sourceinfo = NA_character_, ott_id = ott, node_type = node_type,
      label_status = label_status
    )
  }
  
  nodes <- switch(taxonomy,
                  ncbi = read_ncbi(taxonomy_path),
                  ott = read_ott(taxonomy_path),
                  synthetic = read_newick(taxonomy_path),
                  custom = read_newick(taxonomy_path))
  n <- nrow(nodes)
  if (!n) stop("Scaffold contains no nodes.")
  if (anyNA(nodes$id) || any(!nzchar(nodes$id)) || anyDuplicated(nodes$id)) {
    stop("Scaffold IDs must be nonmissing and unique.")
  }
  if (!"ott_id" %in% names(nodes)) nodes[, ott_id := if (taxonomy == "ott") id else NA_character_]
  if (!"node_type" %in% names(nodes)) nodes[, node_type := "taxon"]
  if (!"label_status" %in% names(nodes)) nodes[, label_status := "native_taxon"]
  if (taxonomy %in% c("synthetic", "custom")) {
    say("Newick labels: ", sum(nodes$label_status == "ott_taxon"),
        " OTT taxa; ", sum(nodes$label_status == "synthetic_mrca"),
        " synthetic MRCAs; ", sum(nodes$label_status == "missing"),
        " missing; ", sum(nodes$label_status == "unrecognized"), " unrecognized.")
  }
  if (length(newick_warnings)) {
    say("ape::read.tree() emitted ", length(newick_warnings),
        " numeric-coercion warnings; verify source labels if mapping is incomplete.")
  }
  
  # For NCBI, the root conventionally points to itself. Convert to NA.
  nodes[parent == id, parent := NA_character_]
  parent_idx <- match(nodes$parent, nodes$id)
  if (any(!is.na(nodes$parent) & is.na(parent_idx))) {
    bad <- nodes$id[!is.na(nodes$parent) & is.na(parent_idx)]
    stop("Parent IDs missing from scaffold; first example: ", bad[1L])
  }
  roots <- which(is.na(parent_idx))
  if (length(roots) != 1L) {
    stop("Expected one rooted connected hierarchy; found ", length(roots), " roots.")
  }
  # Count immediate children and traverse leaves -> root using a queue.
  children_remaining <- tabulate(parent_idx[!is.na(parent_idx)], nbins = n)
  child_count <- children_remaining
  tip_count <- as.numeric(child_count == 0L)
  if (taxonomy == "synthetic" && !is.null(ott_path)) {
    ott_meta <- read_ott(ott_path)
    if (anyDuplicated(ott_meta$id)) stop("Duplicate OTT IDs in rank annotations.")
    idx <- match(nodes$ott_id, ott_meta$id)
    nodes[, rank := ott_meta$rank[idx]]
    nodes[, flags := ott_meta$flags[idx]]
    nodes[, sourceinfo := ott_meta$sourceinfo[idx]]
    matched <- sum(!is.na(idx))
    recognized <- sum(!is.na(nodes$ott_id))
    say("OTT taxonomy matched: ", format(matched, big.mark = ","),
        " / ", format(recognized, big.mark = ","),
        " recognized OTT IDs; species-ranked nodes: ",
        format(sum(nodes$rank == "species", na.rm = TRUE), big.mark = ","))
    if (recognized > 0L && matched == 0L) {
      stop("No extracted OTT IDs matched OTT taxonomy.tsv. Check ott_path and OTT versions.")
    }
  }
  # Rank metadata is absent for a bare custom tree or synthetic tree without OTT.
  has_species_ranks <- taxonomy %in% c("ncbi", "ott") ||
    (taxonomy == "synthetic" && !is.null(ott_path))
  species_count <- if (has_species_ranks) as.numeric(nodes$rank == "species" & !is.na(nodes$rank)) else rep(NA_real_, n)
  if (has_species_ranks && !include_suppressed && taxonomy %in% c("ott", "synthetic")) {
    # OTT flags are comma-separated; 'suppress' marks excluded taxa.
    suppressed <- !is.na(nodes$flags) & grepl("(^|,)\\s*(suppress|hidden|not_otu)(,|$)", nodes$flags)
    species_count[suppressed] <- 0
  }
  depth <- integer(n)
  queue <- integer(n)
  initial <- which(children_remaining == 0L)
  tail <- length(initial)
  queue[seq_len(tail)] <- initial
  head <- 1L
  processed <- 0L
  # First pass: postorder and subtree totals.
  postorder <- integer(n)
  while (head <= tail) {
    i <- queue[head]
    head <- head + 1L
    processed <- processed + 1L
    postorder[processed] <- i
    p <- parent_idx[i]
    if (!is.na(p)) {
      tip_count[p] <- tip_count[p] + tip_count[i]
      if (has_species_ranks) species_count[p] <- species_count[p] + species_count[i]
      children_remaining[p] <- children_remaining[p] - 1L
      if (children_remaining[p] == 0L) {
        tail <- tail + 1L
        queue[tail] <- p
      }
    }
  }
  if (processed != n) stop("Hierarchy contains a cycle or invalid parent relationships.")
  # Second pass: root-first order to compute depth.
  for (i in rev(postorder)) {
    p <- parent_idx[i]
    if (!is.na(p)) depth[i] <- depth[p] + 1L
  }
  nodes[, `:=`(species_count = species_count, tip_count = tip_count,
               depth = depth, child_count = child_count)]
  # Store a dedicated mapping table. Sourceinfo is a comma-separated list of
  # source namespace:identifier tokens; preserve all matches without assuming
  # one-to-one correspondence.
  say("Building identifier cross-references...")
  refs <- data.table::data.table(ott_id = character(), source = character(), source_id = character())
  if (taxonomy == "ott" || (taxonomy == "synthetic" && !is.null(ott_path))) {
    valid <- which(!is.na(nodes$sourceinfo) & nzchar(nodes$sourceinfo) & !is.na(nodes$ott_id))
    if (length(valid)) {
      parts <- strsplit(nodes$sourceinfo[valid], ",", fixed = TRUE)
      lens <- lengths(parts)
      tokens <- trim(unlist(parts, use.names = FALSE))
      oid <- rep(nodes$ott_id[valid], lens)
      ok <- grepl("^[^:]+:.+$", tokens)
      refs <- data.table::data.table(
        ott_id = oid[ok], source = sub(":.*$", "", tokens[ok]),
        source_id = sub("^[^:]+:", "", tokens[ok]))
      refs <- unique(refs)
    }
  }
  # Write to a temporary DB, then rename only after a successful build.
  temp_db <- tempfile(pattern = "taxonsampling_", tmpdir = output_dir, fileext = ".duckdb")
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = temp_db, read_only = FALSE)
  connected <- TRUE
  on.exit({
    if (connected) try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    if (file.exists(temp_db)) unlink(temp_db)
  }, add = TRUE)
  say("Writing DuckDB database...")
  DBI::dbWriteTable(con, "nodes", as.data.frame(nodes), overwrite = TRUE)
  written_cols <- DBI::dbListFields(con, "nodes")
  required_cols <- c("id", "parent", "ott_id", "rank", "node_type", "label_status",
                     "species_count", "tip_count")
  if (!all(required_cols %in% written_cols)) {
    stop("Written nodes table is missing columns: ",
         paste(setdiff(required_cols, written_cols), collapse = ", "))
  }
  DBI::dbWriteTable(con, "crossreferences", as.data.frame(refs), overwrite = TRUE)
  metadata <- data.frame(
    key = c("taxonomy", "source_path", "build_date", "node_count",
            "root_id", "species_count_available", "include_suppressed",
            "species_definition", "tip_definition",
            "ott_mapped_nodes", "missing_node_labels", "unrecognized_node_labels",
            "newick_numeric_warnings"),
    value = c(taxonomy, normalizePath(taxonomy_path), as.character(Sys.Date()),
              as.character(n), nodes$id[roots], as.character(has_species_ranks),
              as.character(include_suppressed),
              if (has_species_ranks) "species-ranked nodes in complete scaffold" else "unavailable: rank annotations missing",
              "terminal nodes in complete scaffold",
              as.character(sum(!is.na(nodes$ott_id))),
              as.character(sum(nodes$label_status == "missing")),
              as.character(sum(nodes$label_status == "unrecognized")),
              as.character(length(newick_warnings))),
    stringsAsFactors = FALSE
  )
  DBI::dbWriteTable(con, "metadata", metadata, overwrite = TRUE)
  DBI::dbExecute(con, "CREATE UNIQUE INDEX nodes_id_idx ON nodes(id)")
  DBI::dbExecute(con, "CREATE INDEX nodes_parent_idx ON nodes(parent)")
  DBI::dbExecute(con, "CREATE INDEX nodes_ott_idx ON nodes(ott_id)")
  DBI::dbExecute(con, "CREATE INDEX refs_ott_idx ON crossreferences(ott_id)")
  DBI::dbExecute(con, "CREATE INDEX refs_source_idx ON crossreferences(source, source_id)")
  DBI::dbDisconnect(con, shutdown = TRUE)
  connected <- FALSE
  if (file.exists(destination)) unlink(destination)
  if (!file.rename(temp_db, destination)) stop("Could not move database to: ", destination)
  say("Saved ", format(n, big.mark = ","), " nodes to ", destination)
  invisible(destination)
}
