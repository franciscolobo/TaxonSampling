# Backbone-independent TaxonSampling engine.
# Requires a taxonsampling object from the revised get_counts().

.ts_draw_one <- function(values, prob = NULL) {
  # sample(x, 1) is unsafe when x is a single numeric-like ID.
  values[[sample.int(length(values), 1L, prob = prob)]]
}

.ts_allocate <- function(m, lower, capacity, weights, method, randomize) {
  stopifnot(length(lower) == length(capacity),
            length(capacity) == length(weights),
            all(lower <= capacity))
  out <- as.integer(lower)
  names(out) <- names(lower)
  remaining <- as.integer(m - sum(out))
  if (remaining < 0L) stop("Required taxa exceed requested sample size.")
  if (is.finite(sum(capacity)) && remaining > sum(capacity - out))
    stop("Requested allocation exceeds available capacity.")
  if (!remaining) return(out)
  if (method == "balanced" && randomize != "yes") {
    # Balanced water filling, with random tie-breaking.
    while (remaining > 0L) {
      eligible <- which(out < capacity)
      if (!length(eligible)) stop("Insufficient capacity.")
      low <- eligible[out[eligible] == min(out[eligible])]
      if (length(low) > remaining) {
        if (length(low) == 1L) chosen <- low else {
          pr <- if (sum(weights[low]) > 0) weights[low] else NULL
          chosen <- low[sample.int(length(low), remaining, prob = pr)]
        }
      } else chosen <- low
      out[chosen] <- out[chosen] + 1L
      remaining <- remaining - length(chosen)
    }
    return(out)
  }
  # Diversity: give each eligible lineage one representative per round.
  # after_first_round: deterministic first coverage, stochastic thereafter.
  first_round <- method == "diversity" && randomize != "yes"
  if (first_round) {
    eligible <- which(out < capacity)
    if (length(eligible) && remaining >= length(eligible)) {
      out[eligible] <- out[eligible] + 1L
      remaining <- remaining - length(eligible)
    }
  }
  while (remaining > 0L) {
    eligible <- which(out < capacity)
    if (!length(eligible)) stop("Insufficient capacity.")
    if (method == "diversity" && randomize == "no") {
      eligible <- eligible[out[eligible] == min(out[eligible])]
    }
    pr <- weights[eligible]
    if (!all(is.finite(pr)) || sum(pr) <= 0) pr <- rep(1, length(eligible))
    chosen <- eligible[sample.int(length(eligible), 1L, prob = pr)]
    out[chosen] <- out[chosen] + 1L
    remaining <- remaining - 1L
  }
  out
}

#' Sample input taxa using a taxonomy backbone
#' @param taxlist Result of get_counts().
#' @param taxon Native backbone node ID (character).
#' @param m Number of observations to select.
#' @param ignoreIDs,requireIDs Native backbone node IDs (not original TaxIDs).
#' @param replacement Allow repeated observations.
#' @param method "diversity" or "balanced".
#' @param randomize "no", "yes", or "after_first_round".
#' @param sampling "agnostic" or "known_species".
#' @param seq_file Optional FASTA file; output sequences selected by seqID.
#' @param out_file Optional FASTA output path (requires seq_file).
#' @param verbose Print messages.
#' @return Updated taxonsampling object with outputNodes, outputIDs,
#'   outputMapping, and optionally outputSeqs.
#' @export
run_TS <- function(taxlist, taxon, m, seq_file = NULL, out_file = NULL,
                   method = c("diversity", "balanced"),
                   randomize = c("no", "yes", "after_first_round"),
                   replacement = FALSE, ignoreIDs = NULL, requireIDs = NULL,
                   sampling = c("agnostic", "known_species"), verbose = TRUE) {
  method <- match.arg(method)
  randomize <- match.arg(randomize)
  sampling <- match.arg(sampling)
  if (method == "balanced" && randomize == "after_first_round")
    stop("after_first_round is supported only for diversity sampling.")
  if (!is.list(taxlist) || !all(c("nodes", "ids_df", "countIDs", "spp_df") %in% names(taxlist)))
    stop("taxlist must be an object returned by get_counts().")
  if (length(taxon) != 1L || is.na(taxon)) stop("taxon must be one node ID.")
  taxon <- as.character(taxon)
  if (length(m) != 1L || is.na(m) || !is.finite(m) || m < 1 || m != floor(m))
    stop("m must be a positive integer.")
  if (!is.logical(replacement) || length(replacement) != 1L || is.na(replacement))
    stop("replacement must be TRUE or FALSE.")
  nodes <- taxlist$nodes
  nodes$id <- as.character(nodes$id)
  nodes$parent <- as.character(nodes$parent)
  if (anyDuplicated(nodes$id)) stop("Duplicate node IDs.")
  if (!taxon %in% nodes$id) stop("Requested root is not an observed backbone node.")
  ids <- taxlist$ids_df
  if (!all(c("taxID", "seqID", "node_id") %in% names(ids)))
    stop("ids_df must have taxID, seqID, and node_id columns.")
  ids$taxID <- as.character(ids$taxID)
  ids$seqID <- as.character(ids$seqID)
  ids$node_id <- as.character(ids$node_id)
  if (anyNA(ids$node_id) || anyNA(ids$taxID) || anyNA(ids$seqID) ||
      any(!nzchar(ids$seqID))) stop("Invalid mapped identifiers or sequence IDs.")
  mode <- taxlist$metadata$count_mode
  if (is.null(mode) || !mode %in% c("input_ids", "backbone_nodes"))
    stop("Unrecognized count_mode in taxlist metadata.")
  parent <- stats::setNames(nodes$parent, nodes$id)
  children <- split(nodes$id[nodes$id != nodes$parent & nodes$parent %in% nodes$id],
                    nodes$parent[nodes$id != nodes$parent & nodes$parent %in% nodes$id])
  descendants <- function(root) {
    result <- root
    frontier <- root
    while (length(frontier)) {
      frontier <- unlist(children[frontier], use.names = FALSE)
      if (!length(frontier)) break
      if (any(frontier %in% result)) stop("Cycle detected in taxonomy.")
      result <- c(result, frontier)
    }
    result
  }
  scope <- descendants(taxon)
  validate_nodes <- function(x, label) {
    x <- unique(as.character(x))
    if (anyNA(x) || any(!x %in% nodes$id))
      stop(label, " contains IDs absent from the observed backbone.")
    if (any(!x %in% scope)) stop(label, " contains IDs outside the selected root.")
    x
  }
  ignored <- validate_nodes(ignoreIDs, "ignoreIDs")
  required <- validate_nodes(requireIDs, "requireIDs")
  # Remove any input observation below an excluded clade.
  excluded <- unique(unlist(lapply(ignored, descendants), use.names = FALSE))
  ids <- ids[ids$node_id %in% scope & !ids$node_id %in% excluded, , drop = FALSE]
  if (!nrow(ids)) stop("No eligible observations remain.")
  if (any(required %in% excluded)) stop("A required clade is excluded.")
  # One row per taxon, or one row per backbone node, is an observation.
  observations <- if (mode == "input_ids") {
    ids[!duplicated(ids$taxID), c("taxID", "node_id"), drop = FALSE]
  } else {
    ids[!duplicated(ids$node_id), c("taxID", "node_id"), drop = FALSE]
  }
  # Requirements specify clade coverage, not individual sequence records.
  # An ancestor requirement is redundant if a required descendant exists.
  if (length(required)) {
    redundant <- vapply(required, function(x)
      any(vapply(setdiff(required, x), function(y) x %in% {
        chain <- y
        while (length(chain) && !is.na(parent[chain]) && parent[chain] != chain) {
          chain <- c(chain, unname(parent[tail(chain, 1L)]))
          if (length(chain) > length(nodes)) stop("Cycle detected.")
        }
        chain
      }, logical(1))), logical(1))
    required <- required[!redundant]
  }
  # Postorder counts and required coverage are computed on the observed tree.
  direct <- tabulate(match(observations$node_id, nodes$id), nbins = nrow(nodes))
  names(direct) <- nodes$id
  bg <- stats::setNames(as.numeric(taxlist$spp_df$species_count),
                        as.character(taxlist$spp_df$taxID))
  if (anyNA(bg)) stop("Background counts contain NA.")
  # Required clades are disjoint after redundant ancestor removal.
  required_below <- function(node) sum(required %in% descendants(node))
  total <- function(node) sum(direct[descendants(node)])
  available <- total(taxon)
  if (!replacement && m > available) {
    warning("Requested m exceeds eligible observations; reducing to ", available, ".")
    m <- available
  }
  if (m < 1L) stop("No observations available after filtering.")
  if (m < length(required)) stop("m is smaller than the number of required clades.")
  for (r in required) if (total(r) < 1) stop("Required clade has no observations: ", r)
  # Recursion selects observed backbone nodes; terminal mapping is separate.
  draw_node <- function(node, n) {
    if (n == 0L) return(character())
    kids <- intersect(children[[node]], scope)
    kids <- kids[vapply(kids, total, numeric(1)) > 0]
    opts <- kids
    if (direct[node] > 0) opts <- c(opts, node)
    if (!length(opts)) stop("No observed taxon below node: ", node)
    cap <- vapply(opts, function(x) if (x == node) direct[node] else total(x), numeric(1))
    if (replacement) cap[] <- Inf
    lower <- vapply(opts, function(x) if (x == node) as.integer(node %in% required)
                    else required_below(x), integer(1))
    # Required ancestor clades are satisfied by any selected descendant.
    if (node %in% required && !any(lower > 0L)) {
      # A requirement at this node needs one sample anywhere below it.
      # n >= 1 is enforced by the parent allocation.
    }
    if (n < sum(lower)) stop("Infeasible required allocation below ", node)
    weights <- if (sampling == "known_species") {
      vapply(opts, function(x) if (x == node) 1 else {
        z <- bg[x]
        if (is.na(z)) 0 else max(0, z)
      }, numeric(1))
    } else rep(1, length(opts))
    names(cap) <- names(lower) <- names(weights) <- opts
    alloc <- .ts_allocate(n, lower, cap, weights, method, randomize)
    result <- character()
    for (x in opts[alloc > 0]) {
      if (x == node) result <- c(result, rep(node, alloc[x]))
      else result <- c(result, draw_node(x, alloc[x]))
    }
    result
  }
  selected_nodes <- draw_node(taxon, as.integer(m))
  if (length(selected_nodes) != m) stop("Internal error: sample size mismatch.")
  # Resolve each sampled node to a unique observation when no replacement.
  remaining <- observations
  output <- vector("list", length(selected_nodes))
  for (i in seq_along(selected_nodes)) {
    node <- selected_nodes[i]
    choices <- which(remaining$node_id == node)
    if (!length(choices)) stop("Internal error: exhausted node ", node)
    j <- choices[sample.int(length(choices), 1L)]
    original <- remaining$taxID[j]
    sequence_choices <- ids$seqID[ids$taxID == original & ids$node_id == node]
    if (!length(sequence_choices)) stop("No sequence ID for ", original)
    seqid <- .ts_draw_one(sequence_choices)
    output[[i]] <- data.frame(sample_index = i, node_id = node,
                              taxID = original, seqID = seqid,
                              stringsAsFactors = FALSE)
    if (!replacement) remaining <- remaining[-j, , drop = FALSE]
  }
  output <- do.call(rbind, output)
  rownames(output) <- NULL
  taxlist$outputNodes <- output$node_id
  taxlist$outputIDs <- output$taxID
  taxlist$outputMapping <- output
  taxlist$ts.params <- list(taxon = taxon, m = m, method = method,
                            randomize = randomize, replacement = replacement,
                            ignoreIDs = ignored, requireIDs = required,
                            sampling = sampling)
  if (!is.null(seq_file)) {
    taxlist <- extract_sequences(taxlist, seq_file, verbose = verbose)
    if (!is.null(out_file)) seqinr::write.fasta(sequences = taxlist$outputSeqs,
                                               names = names(taxlist$outputSeqs),
                                               file.out = out_file)
  } else if (!is.null(out_file)) stop("out_file requires seq_file.")
  if (verbose) message("Selected ", length(taxlist$outputIDs), " observations.")
  taxlist
}
