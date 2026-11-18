#!/usr/bin/env Rscript
# Standalone integration and integrity tests for TaxonSampling cross_identity.duckdb.
# Usage: Rscript test_cross_identity.R taxonomy_data/taxonomy_db [output_directory]
args <- commandArgs(trailingOnly = TRUE)
db_dir <- if (length(args) >= 1L) args[[1L]] else "taxonomy_data/taxonomy_db"
out_dir <- if (length(args) >= 2L) args[[2L]] else file.path(db_dir, "cross_identity_validation")
for (p in c("DBI", "duckdb")) if (!requireNamespace(p, quietly = TRUE)) stop("Install package: ", p)
if (!dir.exists(db_dir)) stop("Missing directory: ", db_dir)
paths <- setNames(file.path(db_dir, paste0(c("ncbi", "ott", "synthetic", "cross_identity"), ".duckdb")), c("ncbi", "ott", "synthetic", "cross_identity"))
if (!file.exists(paths[["cross_identity"]])) stop("Missing cross_identity.duckdb")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
con <- DBI::dbConnect(duckdb::duckdb(), dbdir = normalizePath(paths[["cross_identity"]]), read_only = TRUE)
on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)
query <- function(sql) DBI::dbGetQuery(con, sql)
scalar <- function(sql) as.numeric(query(sql)[[1L]][1L])
quote_sql <- function(x) as.character(DBI::dbQuoteString(con, x))
for (src in c("ncbi", "ott", "synthetic")) {
  if (file.exists(paths[[src]])) {
    DBI::dbExecute(con, paste0("ATTACH ", quote_sql(normalizePath(paths[[src]])), " AS ", src, " (READ_ONLY)"))
  }
}
results <- data.frame(test = character(), passed = logical(), detail = character())
check <- function(name, condition, detail = "") {
  passed <- isTRUE(condition)
  results <<- rbind(results, data.frame(test = name, passed = passed, detail = as.character(detail)))
  cat(if (passed) "PASS " else "FAIL ", name, if (nzchar(detail)) paste0(" — ", detail) else "", "\n", sep = "")
}
safe <- function(name, expr) {
  tryCatch(force(expr), error = function(e) check(name, FALSE, conditionMessage(e)))
}
required <- c("identity_links", "identity_links_status", "transitive_links", "metadata")
tables <- DBI::dbListTables(con)
check("Required cross-identity tables", all(required %in% tables), paste(setdiff(required, tables), collapse = ", "))
if (!all(required %in% tables)) stop("Cannot continue: missing required tables")

safe("No null identifier fields", {
  n <- scalar("SELECT COUNT(*) FROM identity_links_status WHERE source_a IS NULL OR source_b IS NULL OR id_a IS NULL OR id_b IS NULL OR id_a='' OR id_b=''")
  check("No null identifier fields", n == 0, paste("invalid links:", n))
})
safe("No duplicate direct links", {
  n <- scalar("SELECT COUNT(*) FROM (SELECT source_a,source_b,id_a,id_b,COUNT(*) AS n FROM identity_links_status GROUP BY 1,2,3,4 HAVING COUNT(*)>1)")
  check("No duplicate direct links", n == 0, paste("duplicated pairs:", n))
})
safe("Direct mapping status matches degrees", {
  n <- scalar("SELECT COUNT(*) FROM identity_links_status WHERE degree_a < 1 OR degree_b < 1 OR mapping_status <> CASE WHEN degree_a=1 AND degree_b=1 THEN 'unique' ELSE 'ambiguous' END")
  check("Direct mapping status matches degrees", n == 0, paste("inconsistent rows:", n))
})
safe("Stored degrees match recomputed degrees", {
  n <- scalar("WITH a AS (SELECT source_a,source_b,id_a,COUNT(*) AS degree FROM identity_links_status GROUP BY 1,2,3), b AS (SELECT source_a,source_b,id_b,COUNT(*) AS degree FROM identity_links_status GROUP BY 1,2,3) SELECT COUNT(*) FROM identity_links_status l JOIN a USING(source_a,source_b,id_a) JOIN b USING(source_a,source_b,id_b) WHERE l.degree_a <> a.degree OR l.degree_b <> b.degree")
  check("Stored degrees match recomputed degrees", n == 0, paste("inconsistent rows:", n))
})
safe("Unique direct links are reversible", {
  n <- scalar("SELECT COUNT(*) FROM identity_links_status WHERE mapping_status='unique' AND (degree_a<>1 OR degree_b<>1)")
  check("Unique direct links are reversible", n == 0, paste("violations:", n))
})

# Verify links against independently stored native node IDs and source crossreferences.
if (all(file.exists(paths[c("ncbi", "ott")]))) {
  safe("NCBI-OTT links correspond to native nodes and source references", {
    n <- scalar("SELECT COUNT(*) FROM identity_links_status l LEFT JOIN ncbi.nodes n ON n.id=l.id_a LEFT JOIN ott.nodes o ON o.id=l.id_b LEFT JOIN ott.crossreferences x ON x.ott_id=l.id_b AND x.source_id=l.id_a AND lower(x.source)='ncbi' WHERE l.source_a='ncbi' AND l.source_b='ott' AND (n.id IS NULL OR o.id IS NULL OR x.ott_id IS NULL)")
    check("NCBI-OTT links correspond to native nodes and source references", n == 0, paste("invalid links:", n))
  })
}
if (all(file.exists(paths[c("ott", "synthetic")]))) {
  safe("OTT-Synthetic links correspond to native OTT IDs", {
    n <- scalar("SELECT COUNT(*) FROM identity_links_status l LEFT JOIN ott.nodes o ON o.id=l.id_a LEFT JOIN synthetic.nodes s ON s.id=l.id_b WHERE l.source_a='ott' AND l.source_b='synthetic' AND (o.id IS NULL OR s.id IS NULL OR s.ott_id IS DISTINCT FROM l.id_a)")
    check("OTT-Synthetic links correspond to native OTT IDs", n == 0, paste("invalid links:", n))
  })
}

# All transitive rows must be backed by both direct edges; unique paths must be composed of unique edges.
safe("Transitive paths have direct-link evidence", {
  n <- scalar("SELECT COUNT(*) FROM transitive_links t WHERE NOT EXISTS (SELECT 1 FROM identity_links_status a JOIN identity_links_status b ON a.id_b=b.id_a WHERE a.source_a='ncbi' AND a.source_b='ott' AND b.source_a='ott' AND b.source_b='synthetic' AND a.id_a=t.ncbi_id AND a.id_b=t.ott_id AND b.id_b=t.synthetic_id)")
  check("Transitive paths have direct-link evidence", n == 0, paste("invalid paths:", n))
})
safe("Unique transitive paths use unique direct edges", {
  n <- scalar("SELECT COUNT(*) FROM transitive_links t JOIN identity_links_status a ON a.source_a='ncbi' AND a.source_b='ott' AND a.id_a=t.ncbi_id AND a.id_b=t.ott_id JOIN identity_links_status b ON b.source_a='ott' AND b.source_b='synthetic' AND b.id_a=t.ott_id AND b.id_b=t.synthetic_id WHERE t.mapping_status='unique_path' AND (a.mapping_status<>'unique' OR b.mapping_status<>'unique')")
  check("Unique transitive paths use unique direct edges", n == 0, paste("invalid paths:", n))
})
safe("No duplicate transitive paths", {
  n <- scalar("SELECT COUNT(*) FROM (SELECT ncbi_id,ott_id,synthetic_id,COUNT(*) n FROM transitive_links GROUP BY 1,2,3 HAVING COUNT(*)>1)")
  check("No duplicate transitive paths", n == 0, paste("duplicate paths:", n))
})

# Critical: a unique *path* may still have multiple endpoints. This is a diagnostic,
# not necessarily a schema failure; map_taxa(strict) should reject such mappings.
safe("Transitive endpoint ambiguity diagnostic", {
  amb <- query("WITH a AS (SELECT ncbi_id,COUNT(DISTINCT synthetic_id) n FROM transitive_links WHERE mapping_status='unique_path' GROUP BY 1), b AS (SELECT synthetic_id,COUNT(DISTINCT ncbi_id) n FROM transitive_links WHERE mapping_status='unique_path' GROUP BY 1) SELECT COUNT(*) AS ambiguous_unique_paths FROM transitive_links t JOIN a USING(ncbi_id) JOIN b USING(synthetic_id) WHERE t.mapping_status='unique_path' AND (a.n>1 OR b.n>1)")
  write.csv(amb, file.path(out_dir, "transitive_endpoint_ambiguity.csv"), row.names = FALSE)
  check("Transitive endpoint ambiguity diagnostic", TRUE, paste("unique-path rows with nonunique endpoints:", amb$ambiguous_unique_paths[1]))
})

coverage <- query("SELECT source_a,source_b,mapping_status,COUNT(*) AS links, COUNT(DISTINCT id_a) AS distinct_source_a,COUNT(DISTINCT id_b) AS distinct_source_b FROM identity_links_status GROUP BY 1,2,3 ORDER BY 1,2,3")
write.csv(coverage, file.path(out_dir, "direct_mapping_summary.csv"), row.names = FALSE)
transitive <- query("SELECT mapping_status,COUNT(*) AS paths,COUNT(DISTINCT ncbi_id) AS distinct_ncbi,COUNT(DISTINCT ott_id) AS distinct_ott,COUNT(DISTINCT synthetic_id) AS distinct_synthetic FROM transitive_links GROUP BY 1 ORDER BY 1")
write.csv(transitive, file.path(out_dir, "transitive_mapping_summary.csv"), row.names = FALSE)

# Species-only coverage: a crossreference is not necessarily species-to-species.
if (all(file.exists(paths[c("ncbi", "ott")]))) {
  safe("NCBI-OTT species mapping coverage", {
    x <- query("SELECT COUNT(*) AS links,COUNT(DISTINCT l.id_a) AS ncbi_species,COUNT(DISTINCT l.id_b) AS ott_species FROM identity_links_status l JOIN ncbi.nodes n ON n.id=l.id_a JOIN ott.nodes o ON o.id=l.id_b WHERE l.source_a='ncbi' AND l.source_b='ott' AND n.rank='species' AND o.rank='species' AND l.mapping_status='unique'")
    write.csv(x, file.path(out_dir, "ncbi_ott_species_coverage.csv"), row.names = FALSE)
    check("NCBI-OTT species mapping coverage", TRUE, paste("unique species pairs:", x$links[1]))
  })
}
if (all(file.exists(paths[c("ott", "synthetic")]))) {
  safe("OTT-Synthetic species mapping coverage", {
    x <- query("SELECT COUNT(*) AS links,COUNT(DISTINCT l.id_a) AS ott_species,COUNT(DISTINCT l.id_b) AS synthetic_species FROM identity_links_status l JOIN ott.nodes o ON o.id=l.id_a JOIN synthetic.nodes s ON s.id=l.id_b WHERE l.source_a='ott' AND l.source_b='synthetic' AND o.rank='species' AND s.rank='species' AND l.mapping_status='unique'")
    write.csv(x, file.path(out_dir, "ott_synthetic_species_coverage.csv"), row.names = FALSE)
    check("OTT-Synthetic species mapping coverage", TRUE, paste("unique species pairs:", x$links[1]))
  })
}

# Test the public API when the function is available in the R session or source path.
# Supply a second command-line option via TS_CROSS_IDENTITY_R to source the builder.
api_file <- Sys.getenv("TS_CROSS_IDENTITY_R", "")
if (nzchar(api_file)) source(api_file, local = FALSE)
if (exists("map_taxa", mode = "function")) {
  safe("map_taxa strict round-trip for NCBI-OTT", {
    sample_ids <- query("SELECT id_a FROM identity_links_status WHERE source_a='ncbi' AND source_b='ott' AND mapping_status='unique' ORDER BY id_a LIMIT 100")$id_a
    if (!length(sample_ids)) {
      check("map_taxa strict round-trip for NCBI-OTT", FALSE, "No unique NCBI-OTT links")
    } else {
      fw <- map_taxa(sample_ids, "ncbi", "ott", paths[["cross_identity"]], "strict")
      bw <- map_taxa(fw$mapped_id, "ott", "ncbi", paths[["cross_identity"]], "strict")
      joined <- merge(fw, bw, by.x = "mapped_id", by.y = "input_id")
      good <- nrow(joined) == length(sample_ids) && all(joined$input_id == joined$mapped_id.y)
      check("map_taxa strict round-trip for NCBI-OTT", good, paste("tested:", length(sample_ids), "returned:", nrow(joined)))
    }
  })
  safe("map_taxa strict excludes ambiguous direct links", {
    ids <- query("SELECT DISTINCT id_a FROM identity_links_status WHERE source_a='ncbi' AND source_b='ott' AND mapping_status='ambiguous' LIMIT 100")$id_a
    if (!length(ids)) check("map_taxa strict excludes ambiguous direct links", TRUE, "No ambiguous links to test") else {
      out <- map_taxa(ids, "ncbi", "ott", paths[["cross_identity"]], "strict")
      check("map_taxa strict excludes ambiguous direct links", !any(out$input_id %in% ids), paste("tested:", length(ids)))
    }
  })
} else {
  cat("NOTE: map_taxa() API tests skipped; set TS_CROSS_IDENTITY_R to its R source file.\n")
}

write.csv(results, file.path(out_dir, "test_results.csv"), row.names = FALSE)
cat("\nDirect mapping summary:\n"); print(coverage, row.names = FALSE)
cat("\nTransitive mapping summary:\n"); print(transitive, row.names = FALSE)
cat("\nTests: ", sum(results$passed), " passed; ", sum(!results$passed), " failed\n", sep = "")
cat("Reports: ", normalizePath(out_dir), "\n", sep = "")
if (any(!results$passed)) quit(save = "no", status = 1L)
