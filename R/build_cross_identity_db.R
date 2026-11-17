#' Build a separate cross-taxonomy identifier database
#'
#' Preserves all reported links; does not assume mappings are one-to-one.
#' The native taxonomy databases are never modified.
#' @param db_dir Directory containing ncbi.duckdb, ott.duckdb and synthetic.duckdb.
#' @param output_path Destination for the cross-identity database.
#' @param overwrite Replace an existing cross-identity database.
#' @param verbose Display progress.
#' @return Invisibly returns the output path.
#' @export
build_cross_identity_db <- function(db_dir = "taxonomy_db",
                                    output_path = file.path(db_dir, "cross_identity.duckdb"),
                                    overwrite = FALSE, verbose = TRUE) {
  for (pkg in c("DBI", "duckdb")) {
    if (!requireNamespace(pkg, quietly = TRUE)) stop("Required package: ", pkg)
  }
  db_dir <- normalizePath(db_dir, mustWork = TRUE)
  paths <- setNames(file.path(db_dir, paste0(c("ncbi", "ott", "synthetic"), ".duckdb")),
                    c("ncbi", "ott", "synthetic"))
  present <- stats::setNames(file.exists(paths), names(paths))
  if (!present[["ott"]]) stop("ott.duckdb is required for identifier mapping.")
  if (!any(present[c("ncbi", "synthetic")])) {
    stop("At least one of ncbi.duckdb or synthetic.duckdb is required.")
  }
  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  output_path <- normalizePath(output_path, mustWork = FALSE)
  if (file.exists(output_path) && !overwrite) stop("Already exists: ", output_path)
  temp <- tempfile("cross_identity_", tmpdir = dirname(output_path), fileext = ".duckdb")
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir = temp)
  connected <- TRUE
  on.exit({
    if (connected) try(DBI::dbDisconnect(con, shutdown = TRUE), silent = TRUE)
    if (file.exists(temp)) unlink(temp)
  }, add = TRUE)
  run <- function(sql) DBI::dbExecute(con, sql)
  for (src in names(paths)[present]) {
    run(paste0("ATTACH ", as.character(DBI::dbQuoteString(con, normalizePath(paths[[src]]))),
               " AS ", src, " (READ_ONLY)"))
  }
  # A link is evidence of a correspondence, not necessarily a taxonomic-concept identity.
  run("CREATE TABLE identity_links (
    source_a VARCHAR, id_a VARCHAR, source_b VARCHAR, id_b VARCHAR,
    mapping_method VARCHAR, evidence VARCHAR)")
  if (present[["ncbi"]]) {
    run("INSERT INTO identity_links
      SELECT DISTINCT 'ncbi', CAST(n.id AS VARCHAR), 'ott', CAST(o.id AS VARCHAR),
        'ott_sourceinfo', 'ott.crossreferences: ncbi'
      FROM ott.crossreferences x
      JOIN ott.nodes o ON o.id = x.ott_id
      JOIN ncbi.nodes n ON n.id = x.source_id
      WHERE lower(x.source) = 'ncbi'")
  }
  if (present[["synthetic"]]) {
    run("INSERT INTO identity_links
      SELECT DISTINCT 'ott', CAST(o.id AS VARCHAR), 'synthetic', CAST(s.id AS VARCHAR),
        'synthetic_ott_id', 'synthetic.nodes.ott_id'
      FROM synthetic.nodes s JOIN ott.nodes o ON s.ott_id = o.id
      WHERE s.ott_id IS NOT NULL")
  }
  # Degree is computed for each source pair, never globally across namespaces.
  run("CREATE TABLE identity_links_status AS
    WITH a AS (SELECT source_a, source_b, id_a, count(*) AS degree_a
      FROM identity_links GROUP BY 1,2,3),
    b AS (SELECT source_a, source_b, id_b, count(*) AS degree_b
      FROM identity_links GROUP BY 1,2,3)
    SELECT l.*, a.degree_a, b.degree_b,
      CASE WHEN a.degree_a = 1 AND b.degree_b = 1 THEN 'unique'
           ELSE 'ambiguous' END AS mapping_status
    FROM identity_links l
    JOIN a ON l.source_a=a.source_a AND l.source_b=a.source_b AND l.id_a=a.id_a
    JOIN b ON l.source_a=b.source_a AND l.source_b=b.source_b AND l.id_b=b.id_b")
  if (present[["ncbi"]] && present[["synthetic"]]) {
    run("CREATE TABLE transitive_links AS
      SELECT DISTINCT a.id_a AS ncbi_id, a.id_b AS ott_id, b.id_b AS synthetic_id,
        CASE WHEN a.mapping_status='unique' AND b.mapping_status='unique'
          THEN 'unique_path' ELSE 'ambiguous_path' END AS mapping_status
      FROM identity_links_status a JOIN identity_links_status b ON a.id_b=b.id_a
      WHERE a.source_a='ncbi' AND a.source_b='ott'
        AND b.source_a='ott' AND b.source_b='synthetic'")
  } else {
    run("CREATE TABLE transitive_links (ncbi_id VARCHAR, ott_id VARCHAR,
      synthetic_id VARCHAR, mapping_status VARCHAR)")
  }
  run("CREATE TABLE metadata (key VARCHAR, value VARCHAR)")
  meta <- data.frame(key=c("schema_version", "build_date", "ncbi_available", "ott_available",
                            "synthetic_available", "identity_interpretation"),
                     value=c("1", as.character(Sys.Date()), as.character(present[["ncbi"]]),
                             "TRUE", as.character(present[["synthetic"]]),
                             "cross-referenced identifiers; not guaranteed equivalent taxonomic concepts"))
  DBI::dbWriteTable(con, "metadata", meta, append = TRUE)
  run("CREATE INDEX identity_a_idx ON identity_links_status(source_a, source_b, id_a)")
  run("CREATE INDEX identity_b_idx ON identity_links_status(source_a, source_b, id_b)")
  run("CREATE INDEX transitive_ncbi_idx ON transitive_links(ncbi_id)")
  counts <- DBI::dbGetQuery(con, "SELECT source_a, source_b, mapping_status, count(*) AS links
                                    FROM identity_links_status GROUP BY 1,2,3 ORDER BY 1,2,3")
  DBI::dbDisconnect(con, shutdown = TRUE)
  connected <- FALSE
  if (file.exists(output_path)) unlink(output_path)
  if (!file.rename(temp, output_path)) stop("Could not install ", output_path)
  if (verbose) { message("Created ", output_path); print(counts) }
  invisible(output_path)
}

#' Map native identifiers using cross-identity evidence
#' @param ids Character vector of native identifiers.
#' @param from Source taxonomy: ncbi, ott, synthetic.
#' @param to Destination taxonomy.
#' @param cross_identity_db Cross-identity DuckDB path.
#' @param mode strict (unambiguous paths only) or all.
#' @return A data.frame of candidate mappings and mapping status.
#' @export
map_taxa <- function(ids, from, to,
                     cross_identity_db = "taxonomy_db/cross_identity.duckdb",
                     mode = c("strict", "all")) {
  mode <- match.arg(mode)
  allowed <- c("ncbi", "ott", "synthetic")
  if (length(from)!=1L || length(to)!=1L || !from %in% allowed || !to %in% allowed)
    stop("from and to must be one of: ", paste(allowed, collapse=", "))
  if (!file.exists(cross_identity_db)) stop("Missing cross-identity database")
  ids <- unique(as.character(ids))
  if (!length(ids)) return(data.frame(input_id=character(), mapped_id=character(),
                                      mapping_status=character(), stringsAsFactors=FALSE))
  con <- DBI::dbConnect(duckdb::duckdb(), dbdir=cross_identity_db, read_only=TRUE)
  on.exit(DBI::dbDisconnect(con, shutdown=TRUE), add=TRUE)
  id_values <- paste(vapply(ids, function(x) as.character(DBI::dbQuoteString(con, x)), character(1)), collapse=", ")
  if (from == to) {
    return(data.frame(input_id=ids, mapped_id=ids, mapping_status="same_source"))
  }
  if (setequal(c(from,to), c("ncbi","synthetic"))) {
    a <- if (from=="ncbi") "ncbi_id" else "synthetic_id"
    b <- if (to=="ncbi") "ncbi_id" else "synthetic_id"
    sql <- sprintf("SELECT DISTINCT t.%s AS input_id, t.%s AS mapped_id,
      t.mapping_status FROM transitive_links t WHERE t.%s IN (%s)", a, b, a, id_values)
    out <- DBI::dbGetQuery(con, sql)
    # A unique path is not necessarily a unique endpoint: enforce endpoint degrees.
    if (nrow(out)) {
      out$mapping_status <- ifelse(out$mapping_status=="unique_path", "unique", "ambiguous")
      forward <- table(out$input_id)
      reverse <- table(out$mapped_id)
      out$mapping_status[forward[out$input_id] > 1 | reverse[out$mapped_id] > 1] <- "ambiguous"
    }
  } else {
    forward <- (from=="ncbi" && to=="ott") || (from=="ott" && to=="synthetic")
    a <- if (forward) "id_a" else "id_b"
    b <- if (forward) "id_b" else "id_a"
    sql <- sprintf("SELECT DISTINCT l.%s AS input_id, l.%s AS mapped_id,
      l.mapping_status FROM identity_links_status l
      WHERE l.%s IN (%s) AND l.source_a=%s AND l.source_b=%s", a, b, a, id_values,
      as.character(DBI::dbQuoteString(con, if (forward) from else to)),
      as.character(DBI::dbQuoteString(con, if (forward) to else from)))
    out <- DBI::dbGetQuery(con, sql)
  }
  if (mode=="strict") out <- out[out$mapping_status=="unique", , drop=FALSE]
  out
}
