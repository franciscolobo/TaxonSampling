
#' Count input taxa on a precomputed TaxonSampling scaffold
#'
#' Each distinct input identifier contributes one observation by default.
#' Alternatively, count each distinct mapped backbone node once.
#' Background counts are read from the immutable DuckDB scaffold.
#'
#' @param taxonomy_db Path to a DuckDB created by build_taxonomy_db().
#' @param ids_file Headerless TSV (taxon ID, sequence ID).
#' @param ids_df Data frame with taxon ID and sequence ID.
#' @param id_type "ncbi", "ott", or "native".
#' @param count_mode "input_ids" or "backbone_nodes".
#' @param background_unit "species" or "tips".
#' @param duplicate_ids "first" or "error".
#' @param missing_ids "warn", "error", or "drop".
#' @param spp_file Optional TSV overriding background counts.
#' @param spp_df Optional data frame overriding background counts.
#' @param start_from_species Deprecated compatibility parameter.
#' @param verbose Print progress messages.
#'
#' @return A taxonsampling list containing nodes, ids_df, countIDs,
#'   spp_df, mapping, unmapped_ids, and metadata.
#'
#' @export
get_counts <- function(
    taxonomy_db,
    ids_file = NULL,
    ids_df = NULL,
    id_type = c("ncbi", "ott", "native"),
    count_mode = c("input_ids", "backbone_nodes"),
    background_unit = c("species", "tips"),
    duplicate_ids = c("first", "error"),
    missing_ids = c("warn", "error", "drop"),
    spp_file = NULL,
    spp_df = NULL,
    start_from_species = FALSE,
    verbose = TRUE
) {
  
  id_type <- match.arg(id_type)
  count_mode <- match.arg(count_mode)
  background_unit <- match.arg(background_unit)
  duplicate_ids <- match.arg(duplicate_ids)
  missing_ids <- match.arg(missing_ids)
  
  if (!requireNamespace("DBI", quietly = TRUE) ||
      !requireNamespace("duckdb", quietly = TRUE)) {
    stop("Packages DBI and duckdb are required.")
  }
  
  if (!is.character(taxonomy_db) ||
      length(taxonomy_db) != 1L ||
      is.na(taxonomy_db) ||
      !file.exists(taxonomy_db)) {
    stop("taxonomy_db must be an existing DuckDB file.")
  }
  
  if (!is.logical(verbose) ||
      length(verbose) != 1L ||
      is.na(verbose)) {
    stop("verbose must be TRUE or FALSE.")
  }
  
  if (!is.logical(start_from_species) ||
      length(start_from_species) != 1L ||
      is.na(start_from_species)) {
    stop("start_from_species must be TRUE or FALSE.")
  }
  
  if (start_from_species) {
    warning(
      "start_from_species is deprecated; ",
      "background counts are precomputed.",
      call. = FALSE
    )
  }
  
  if (is.null(ids_file) && is.null(ids_df)) {
    stop("Provide ids_file or ids_df.")
  }
  
  say <- function(...) {
    if (verbose) message(...)
  }
  
  # ---------------------------------------------------------------
  # 1. Input loading
  # ---------------------------------------------------------------
  
  read_pairs <- function(file, df, names_out) {
    
    if (!is.null(file)) {
      
      if (length(file) != 1L ||
          is.na(file) ||
          !file.exists(file)) {
        stop("File not found: ", file)
      }
      
      df <- utils::read.delim(
        file,
        header = FALSE,
        sep = "\t",
        quote = "",
        comment.char = "",
        colClasses = "character",
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }
    
    if (!is.data.frame(df) ||
        ncol(df) < 2L ||
        nrow(df) == 0L) {
      stop(
        "Expected a nonempty two-column data frame for ",
        names_out[1L]
      )
    }
    
    df <- as.data.frame(
      df[, 1:2, drop = FALSE],
      stringsAsFactors = FALSE
    )
    
    names(df) <- names_out
    
    df[[1L]] <- trimws(as.character(df[[1L]]))
    df[[2L]] <- as.character(df[[2L]])
    
    df
  }
  
  input <- read_pairs(
    ids_file,
    ids_df,
    c("taxID", "seqID")
  )
  
  invalid <- is.na(input$taxID) |
    !nzchar(input$taxID)
  
  if (any(invalid)) {
    
    warning(
      sum(invalid),
      " rows with empty taxon IDs were discarded.",
      call. = FALSE
    )
    
    input <- input[!invalid, , drop = FALSE]
  }
  
  if (!nrow(input)) {
    stop("No valid input identifiers remain.")
  }
  
  if (id_type == "ncbi") {
    input$taxID <- sub("^ncbi:", "", input$taxID)
  }
  
  if (id_type == "ott") {
    input$taxID <- sub("^ott:", "", input$taxID)
  }
  
  duplicated_input <- duplicated(input$taxID)
  
  if (any(duplicated_input)) {
    
    if (duplicate_ids == "error") {
      stop("Duplicated input taxonomic identifiers.")
    }
    
    warning(
      sum(duplicated_input),
      " duplicated input taxon IDs: ",
      "keeping first occurrence.",
      call. = FALSE
    )
    
    input <- input[
      !duplicated_input,
      ,
      drop = FALSE
    ]
  }
  
  rownames(input) <- NULL
  
  say("Distinct input taxa: ", nrow(input))
  
  # ---------------------------------------------------------------
  # 2. Open reference database
  # ---------------------------------------------------------------
  
  con <- DBI::dbConnect(
    duckdb::duckdb(),
    dbdir = taxonomy_db,
    read_only = TRUE
  )
  
  on.exit(
    DBI::dbDisconnect(
      con,
      shutdown = TRUE
    ),
    add = TRUE
  )
  
  tables <- DBI::dbListTables(con)
  
  if (!all(c("nodes", "metadata") %in% tables)) {
    stop(
      "Invalid taxonomy database: ",
      "nodes or metadata table missing."
    )
  }
  
  md <- DBI::dbGetQuery(
    con,
    "SELECT key, value FROM metadata"
  )
  
  db_metadata <- stats::setNames(
    as.character(md$value),
    md$key
  )
  
  taxonomy <- unname(db_metadata["taxonomy"])
  
  if (length(taxonomy) != 1L ||
      is.na(taxonomy)) {
    taxonomy <- "unknown"
  }
  
  say("Taxonomy: ", taxonomy)
  
  # ---------------------------------------------------------------
  # 3. Load scaffold
  # ---------------------------------------------------------------
  
  say("Loading taxonomy nodes...")
  
  nodes_all <- DBI::dbGetQuery(
    con,
    paste(
      "SELECT id, parent, name, rank, flags,",
      "sourceinfo, ott_id, node_type,",
      "species_count, tip_count, depth, child_count",
      "FROM nodes"
    )
  )
  
  nodes_all$id <- as.character(nodes_all$id)
  nodes_all$parent <- as.character(nodes_all$parent)
  nodes_all$ott_id <- as.character(nodes_all$ott_id)
  
  n <- nrow(nodes_all)
  
  if (!n ||
      anyNA(nodes_all$id) ||
      anyDuplicated(nodes_all$id)) {
    stop(
      "Scaffold node identifiers are ",
      "missing or duplicated."
    )
  }
  
  if (anyNA(nodes_all$depth)) {
    stop("Scaffold contains missing node depths.")
  }
  
  parent_idx <- match(
    nodes_all$parent,
    nodes_all$id
  )
  
  root_idx <- which(
    is.na(nodes_all$parent) |
      !nzchar(nodes_all$parent) |
      nodes_all$parent == nodes_all$id
  )
  
  if (length(root_idx) != 1L) {
    stop("Expected exactly one root node.")
  }
  
  parent_idx[root_idx] <- NA_integer_
  
  if (anyNA(parent_idx[-root_idx])) {
    stop(
      "Some non-root parents are absent ",
      "from scaffold."
    )
  }
  
  if (any(
    nodes_all$depth[-root_idx] <=
    nodes_all$depth[parent_idx[-root_idx]]
  )) {
    stop(
      "Invalid scaffold: child depths must ",
      "exceed parent depths."
    )
  }
  
  # ---------------------------------------------------------------
  # 4. Map input identifiers
  # ---------------------------------------------------------------
  
  say("Mapping input identifiers...")
  
  mapping <- data.frame(
    taxID = input$taxID,
    seqID = input$seqID,
    node_id = NA_character_,
    status = "unmapped",
    stringsAsFactors = FALSE
  )
  
  if (id_type == "native") {
    
    idx <- match(
      mapping$taxID,
      nodes_all$id
    )
    
    mapping$node_id <- nodes_all$id[idx]
    
  } else if (id_type == "ott") {
    
    candidate <- which(
      !is.na(nodes_all$ott_id) &
        nzchar(nodes_all$ott_id)
    )
    
    hits <- candidate[
      match(
        mapping$taxID,
        nodes_all$ott_id[candidate]
      )
    ]
    
    mapping$node_id <- nodes_all$id[hits]
    
    duplicated_ott <- unique(
      nodes_all$ott_id[candidate][
        duplicated(nodes_all$ott_id[candidate])
      ]
    )
    
    mapping$status[
      mapping$taxID %in% duplicated_ott
    ] <- "ambiguous"
    
    mapping$node_id[
      mapping$status == "ambiguous"
    ] <- NA_character_
    
  } else if (taxonomy == "ncbi") {
    
    idx <- match(
      paste0("ncbi:", mapping$taxID),
      nodes_all$id
    )
    
    missing <- is.na(idx)
    
    idx[missing] <- match(
      mapping$taxID[missing],
      nodes_all$id
    )
    
    mapping$node_id <- nodes_all$id[idx]
    
  } else {
    
    if (!"crossreferences" %in% tables) {
      stop(
        "NCBI mapping requires a ",
        "crossreferences table."
      )
    }
    
    requested <- unique(mapping$taxID)
    
    chunks <- split(
      requested,
      ceiling(seq_along(requested) / 1000L)
    )
    
    parts <- lapply(chunks, function(ids) {
      
      quoted <- as.character(
        DBI::dbQuoteString(con, ids)
      )
      
      DBI::dbGetQuery(
        con,
        paste0(
          "SELECT source_id, ott_id ",
          "FROM crossreferences ",
          "WHERE source = 'ncbi' ",
          "AND source_id IN (",
          paste(quoted, collapse = ","),
          ")"
        )
      )
    })
    
    xrefs <- do.call(rbind, parts)
    
    if (nrow(xrefs)) {
      
      xrefs$source_id <- as.character(
        xrefs$source_id
      )
      
      xrefs$ott_id <- as.character(
        xrefs$ott_id
      )
      
      candidates <- which(
        !is.na(nodes_all$ott_id) &
          nzchar(nodes_all$ott_id)
      )
      
      matched <- match(
        xrefs$ott_id,
        nodes_all$ott_id[candidates]
      )
      
      xrefs$node_idx <- candidates[matched]
      
      xrefs <- unique(
        xrefs[
          !is.na(xrefs$node_idx),
          c("source_id", "node_idx"),
          drop = FALSE
        ]
      )
      
      if (nrow(xrefs)) {
        
        target_counts <- table(
          xrefs$source_id
        )
        
        ambiguous <- names(target_counts)[
          target_counts > 1L
        ]
        
        resolved <- xrefs[
          !(xrefs$source_id %in% ambiguous),
          ,
          drop = FALSE
        ]
        
        idx <- match(
          mapping$taxID,
          resolved$source_id
        )
        
        mapping$node_id <- nodes_all$id[
          resolved$node_idx[idx]
        ]
        
        mapping$status[
          mapping$taxID %in% ambiguous
        ] <- "ambiguous"
      }
    }
  }
  
  mapping$status[
    !is.na(mapping$node_id)
  ] <- "mapped"
  
  unmapped_ids <- mapping[
    mapping$status != "mapped",
    ,
    drop = FALSE
  ]
  
  if (nrow(unmapped_ids)) {
    
    msg <- paste0(
      nrow(unmapped_ids),
      " input identifiers could not be mapped uniquely."
    )
    
    if (missing_ids == "error") {
      stop(msg)
    }
    
    if (missing_ids == "warn") {
      warning(msg, call. = FALSE)
    }
  }
  
  mapped <- mapping[
    mapping$status == "mapped",
    ,
    drop = FALSE
  ]
  
  if (!nrow(mapped)) {
    stop("No input identifiers mapped to scaffold.")
  }
  
  say("Mapped input taxa: ", nrow(mapped))
  say("Unmapped/ambiguous taxa: ", nrow(unmapped_ids))
  
  # ---------------------------------------------------------------
  # 5. Initialize counts
  # ---------------------------------------------------------------
  
  mapped_idx <- match(
    mapped$node_id,
    nodes_all$id
  )
  
  count_idx <- if (count_mode == "backbone_nodes") {
    unique(mapped_idx)
  } else {
    mapped_idx
  }
  
  observed <- as.numeric(
    tabulate(
      count_idx,
      nbins = n
    )
  )
  
  say(
    "Counting mode: ",
    count_mode,
    " (",
    length(count_idx),
    " contributions)"
  )
  
  # ---------------------------------------------------------------
  # 6. Build induced ancestor subgraph
  # ---------------------------------------------------------------
  
  active <- rep(FALSE, n)
  
  frontier <- unique(mapped_idx)
  
  while (length(frontier)) {
    
    frontier <- frontier[!active[frontier]]
    
    if (!length(frontier)) {
      break
    }
    
    active[frontier] <- TRUE
    
    frontier <- unique(
      parent_idx[frontier]
    )
    
    frontier <- frontier[
      !is.na(frontier)
    ]
  }
  
  active_idx <- which(active)
  
  say(
    "Nodes in induced taxonomy: ",
    length(active_idx)
  )
  
  # ---------------------------------------------------------------
  # 7. Propagate counts to ancestors
  # ---------------------------------------------------------------
  
  say("Propagating observed counts...")
  
  depths <- sort(
    unique(nodes_all$depth[active_idx]),
    decreasing = TRUE
  )
  
  for (d in depths) {
    
    current <- active_idx[
      nodes_all$depth[active_idx] == d
    ]
    
    parents <- parent_idx[current]
    
    valid <- !is.na(parents)
    
    if (!any(valid)) {
      next
    }
    
    totals <- rowsum(
      observed[current[valid]],
      group = parents[valid],
      reorder = FALSE
    )
    
    targets <- as.integer(
      rownames(totals)
    )
    
    observed[targets] <- observed[targets] +
      as.numeric(totals[, 1L])
  }
  
  # ---------------------------------------------------------------
  # 8. Extract nodes with observed counts
  # ---------------------------------------------------------------
  
  keep <- active_idx[
    observed[active_idx] > 0
  ]
  
  nodes <- nodes_all[
    keep,
    ,
    drop = FALSE
  ]
  
  rownames(nodes) <- NULL
  
  countIDs <- stats::setNames(
    observed[keep],
    nodes$id
  )
  
  # ---------------------------------------------------------------
  # 9. Background counts
  # ---------------------------------------------------------------
  
  background_column <- if (
    background_unit == "species"
  ) {
    "species_count"
  } else {
    "tip_count"
  }
  
  background <- data.frame(
    taxID = nodes$id,
    species_count = as.numeric(
      nodes[[background_column]]
    ),
    stringsAsFactors = FALSE
  )
  
  # ---------------------------------------------------------------
  # 10. Optional background overrides
  # ---------------------------------------------------------------
  
  if (!is.null(spp_file) || !is.null(spp_df)) {
    
    override <- read_pairs(
      spp_file,
      spp_df,
      c("taxID", "species_count")
    )
    
    override$species_count <- suppressWarnings(
      as.numeric(override$species_count)
    )
    
    if (anyNA(override$species_count) ||
        any(!is.finite(override$species_count)) ||
        any(override$species_count < 0) ||
        anyDuplicated(override$taxID)) {
      stop(
        "Invalid or duplicated background overrides."
      )
    }
    
    if (taxonomy == "ncbi") {
      
      raw <- !grepl(
        ":",
        override$taxID,
        fixed = TRUE
      )
      
      override$taxID[raw] <- paste0(
        "ncbi:",
        override$taxID[raw]
      )
    }
    
    idx <- match(
      override$taxID,
      background$taxID
    )
    
    found <- !is.na(idx)
    
    background$species_count[idx[found]] <-
      override$species_count[found]
    
    if (any(!found)) {
      warning(
        sum(!found),
        " overrides did not match retained nodes.",
        call. = FALSE
      )
    }
  }
  
  # ---------------------------------------------------------------
  # 11. Assemble output
  # ---------------------------------------------------------------
  
  result <- list(
    
    nodes = nodes,
    
    ids_df = mapped[
      ,
      c("taxID", "seqID", "node_id"),
      drop = FALSE
    ],
    
    countIDs = countIDs,
    
    spp_df = background,
    
    mapping = mapping,
    
    unmapped_ids = unmapped_ids,
    
    metadata = list(
      
      taxonomy = taxonomy,
      
      taxonomy_db = normalizePath(
        taxonomy_db
      ),
      
      id_type = id_type,
      
      count_mode = count_mode,
      
      background_unit = background_unit,
      
      input_taxa = nrow(input),
      
      mapped_taxa = nrow(mapped),
      
      unique_mapped_nodes = length(
        unique(mapped_idx)
      ),
      
      counted_contributions = length(
        count_idx
      ),
      
      unmapped_taxa = nrow(unmapped_ids),
      
      retained_nodes = nrow(nodes),
      
      root_id = nodes_all$id[root_idx],
      
      root_observed_count = observed[root_idx],
      
      root_background_count =
        background$species_count[
          match(
            nodes_all$id[root_idx],
            background$taxID
          )
        ],
      
      database_metadata = db_metadata
    )
  )
  
  class(result) <- c(
    "taxonsampling",
    "list"
  )
  
  say("Completed get_counts().")
  
  result
}
