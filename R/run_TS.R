# ============================================================================
# TaxonSampling: run_TS.R
#
# High-level flow:
#   1. Validate options and input tables.
#   2. Build an observed-tree adjacency representation.
#   3. Restrict the tree to the requested root and apply clade constraints.
#   4. Count distinct eligible native nodes (capacity) and read background
#      species counts (allocation weights).
#   5. Recursively allocate the requested number of native IDs to branches.
#   6. Map each selected native ID back to a source TaxID and sequence ID.
#   7. Optionally extract/write FASTA sequences and return enriched taxlist.
#
# IMPORTANT: .ts_allocate() delegates the mathematical allocation rules to
# .ts_legacy_fixed_allocate() in R/aaa_legacy_allocation.R.
# This file controls recursion and construction of the allocator inputs.
# ============================================================================

# Backbone-independent TaxonSampling engine.
# Requires a taxonsampling object from the revised get_counts().

# Select one element by its positional index. This avoids base::sample() interpreting a single
# numeric-looking TaxID as an integer range.
.ts_draw_one <- function(values, prob = NULL) {
  # The issue arises when a vector contains exactly one numeric-like identifier.
  # sample(x, 1) is unsafe when x is a single numeric-like ID.
  # sample.int() returns an index; [[ ]] retrieves the original value without coercing its
  # identity.
  values[[sample.int(length(values), 1L, prob = prob)]]
}

# Adapter to the corrected ORIGINAL allocation strategies.
# Adapter used by recursive sampling. m is the TOTAL allocation requested for this node; lower
# contains mandatory allocations already reserved for child branches.
.ts_allocate <- function(m, lower, capacity, weights, method, randomize,
                         sampling = "known_species") {
  # All child-aligned vectors must have the same number of entries.
  if (length(lower) != length(capacity) || length(lower) != length(weights))
    stop("Allocation vectors have different lengths.")
  # Mandatory allocations cannot exceed the total number requested.
  if (m < sum(lower)) stop("Required taxa exceed requested sample size.")
  # The historical .ts_allocate() interface accepts unnamed positional vectors.
  # The legacy allocator requires matching branch names; supply positional IDs
  # only when all three vectors are unnamed, and reject mixed naming.
  # Capture the three naming schemes so the legacy allocator can match branch IDs reliably.
  nms <- list(names(lower), names(capacity), names(weights))
  # For purely positional inputs, generate stable temporary names in the same order for all
  # vectors.
  if (all(vapply(nms, is.null, logical(1)))) {
    branch_ids <- paste0("branch_", seq_along(capacity))
    names(lower) <- names(capacity) <- names(weights) <- branch_ids
    # Mixed named and unnamed inputs are ambiguous and are rejected.
  } else if (any(vapply(nms, is.null, logical(1)))) {
    stop("Allocation vectors must either all be named or all be unnamed.")
  }
  # Delegate the actual allocation formula to aaa_legacy_allocation.R. This file does not itself
  # implement the ten allocation policies.
  .ts_legacy_fixed_allocate(
    # Subtract the mandatory lower bounds because the legacy helper expects only the REMAINING
    # budget; pass those bounds separately as m_i.
    m = m - sum(lower), m_i = lower,
    # Capacity is the maximum count per branch; weights are background diversity counts (or
    # uniform values).
    childrenCount = capacity, childrenCountSpp = weights,
    method = method, randomize = randomize, sampling = sampling)
}

# Public interface. The documentation immediately below is processed by roxygen2.
#' Sample input taxa using a taxonomy backbone
#' @param taxlist Result of get_counts().
#' @param taxon Native backbone node ID (character).
#' @param m Number of distinct taxonomic nodes to select (without replacement).
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
# Main entry point. The argument defaults are ordered vectors because match.arg() uses the first
# option when unspecified.
run_TS <- function(taxlist, taxon, m, seq_file = NULL, out_file = NULL,
                   method = c("diversity", "balanced"),
                   randomize = c("no", "yes", "after_first_round"),
                   replacement = FALSE, ignoreIDs = NULL, requireIDs = NULL,
                   sampling = c("agnostic", "known_species"), verbose = TRUE) {
  # Resolve the user-provided method to exactly one supported string.
  method <- match.arg(method)
  # Resolve randomization mode to exactly one supported string.
  randomize <- match.arg(randomize)
  # Resolve whether allocation weights are uniform or informed by background species diversity.
  sampling <- match.arg(sampling)
  # Balanced allocation has no implemented after-first-round policy; fail early.
  if (method == "balanced" && randomize == "after_first_round")
    stop("after_first_round is supported only for diversity sampling.")
  # Require the structural tables produced by get_counts(), rather than accepting an arbitrary
  # list.
  if (!is.list(taxlist) || !all(c("nodes", "ids_df", "countIDs", "spp_df") %in% names(taxlist)))
    stop("taxlist must be an object returned by get_counts().")
  # The requested root must be a single nonmissing native backbone ID.
  if (length(taxon) != 1L || is.na(taxon)) stop("taxon must be one node ID.")
  # Normalize IDs to character to prevent numeric conversions and preserve opaque IDs such as
  # ott:123.
  taxon <- as.character(taxon)
  # Sample size must be one finite, positive, integer-valued number.
  if (length(m) != 1L || is.na(m) || !is.finite(m) || m < 1 || m != floor(m))
    stop("m must be a positive integer.")
  # Replacement is a single nonmissing Boolean.
  if (!is.logical(replacement) || length(replacement) != 1L || is.na(replacement))
    stop("replacement must be TRUE or FALSE.")
  # Work with the observed taxonomy nodes stored in taxlist.
  nodes <- taxlist$nodes
  # Normalize both ends of every parent-child relationship to character IDs.
  nodes$id <- as.character(nodes$id)
  nodes$parent <- as.character(nodes$parent)
  # Each node must have exactly one record in the backbone table.
  if (anyDuplicated(nodes$id)) stop("Duplicate node IDs.")
  # Only roots represented in this observed backbone are allowed.
  if (!taxon %in% nodes$id) stop("Requested root is not an observed backbone node.")
  # Input mappings connect source TaxIDs and sequence IDs to native backbone nodes.
  ids <- taxlist$ids_df
  # The three columns are required for later mapping of selected nodes back to source data.
  if (!all(c("taxID", "seqID", "node_id") %in% names(ids)))
    stop("ids_df must have taxID, seqID, and node_id columns.")
  # Normalize all identifiers before joins, filtering, and equality checks.
  ids$taxID <- as.character(ids$taxID)
  ids$seqID <- as.character(ids$seqID)
  ids$node_id <- as.character(ids$node_id)
  # Missing mapping identifiers or empty sequence IDs would prevent reliable output
  # reconstruction.
  if (anyNA(ids$node_id) || anyNA(ids$taxID) || anyNA(ids$seqID) ||
      any(!nzchar(ids$seqID))) stop("Invalid mapped identifiers or sequence IDs.")
  # The metadata records how get_counts() interpreted the original input.
  mode <- taxlist$metadata$count_mode
  # Reject unknown modes instead of silently applying incorrect duplicate-handling assumptions.
  if (is.null(mode) || !mode %in% c("input_ids", "backbone_nodes"))
    stop("Unrecognized count_mode in taxlist metadata.")
  # Create a named lookup: parent[[node_id]] yields the parent node ID.
  parent <- stats::setNames(nodes$parent, nodes$id)
  # Build adjacency lists by grouping each non-self-parented child under its parent; ignore
  # parent IDs absent from the observed table.
  children <- split(nodes$id[nodes$id != nodes$parent & nodes$parent %in% nodes$id],
                    nodes$parent[nodes$id != nodes$parent & nodes$parent %in% nodes$id])
  # Helper returns a node together with every descendant reachable through the observed
  # adjacency lists.
  descendants <- function(root) {
    # Initialize the visited set and the current breadth-first frontier with the requested root.
    result <- root
    frontier <- root
    # Continue breadth-first traversal until no new descendants remain.
    while (length(frontier)) {
      # Look up children of every frontier node, flattening the resulting list.
      frontier <- unlist(children[frontier], use.names = FALSE)
      # Stop traversal once a level contains no additional nodes.
      if (!length(frontier)) break
      # Detect repeated visits, which indicate a cycle or otherwise malformed hierarchy.
      if (any(frontier %in% result)) stop("Cycle detected in taxonomy.")
      # Accumulate all nodes discovered at the current level.
      result <- c(result, frontier)
    }
    result
  }
  # Restrict the entire run to the requested root and its descendants.
  scope <- descendants(taxon)
  # Validate user-supplied required and ignored native IDs within this root.
  validate_nodes <- function(x, label) {
    # Convert IDs to character and discard duplicates so repeated constraints have no extra
    # effect.
    x <- unique(as.character(x))
    # Reject missing IDs or IDs absent from the observed backbone.
    if (anyNA(x) || any(!x %in% nodes$id))
      stop(label, " contains IDs absent from the observed backbone.")
    # Constraints cannot reference nodes outside the subtree selected by taxon.
    if (any(!x %in% scope)) stop(label, " contains IDs outside the selected root.")
    x
  }
  # Resolve excluded and required node sets independently.
  ignored <- validate_nodes(ignoreIDs, "ignoreIDs")
  required <- validate_nodes(requireIDs, "requireIDs")
  # Ignoring a clade excludes every descendant, not merely its named root.
  # Remove any input observation below an excluded clade.
  # Expand all ignored roots and deduplicate their descendant nodes.
  excluded <- unique(unlist(lapply(ignored, descendants), use.names = FALSE))
  # Retain mappings that lie inside the selected root but outside all excluded clades.
  ids <- ids[ids$node_id %in% scope & !ids$node_id %in% excluded, , drop = FALSE]
  # Sampling cannot proceed without at least one eligible input mapping.
  if (!nrow(ids)) stop("No eligible observations remain.")
  # A required clade cannot overlap an excluded clade.
  if (any(required %in% excluded)) stop("A required clade is excluded.")
  # IMPORTANT MODEL CHOICE: one distinct native node is one sampling unit when replacement is
  # FALSE.
  # The sampling unit is a distinct native taxonomy node, not a sequence or
  # an input-table row. Multiple input TaxIDs mapping to the same node cannot
  # increase that node's capacity. Select one representative mapping only
  # after the node itself has been sampled.
  # In input-ID mode, verify that a source TaxID does not ambiguously map to multiple native
  # nodes.
  if (mode == "input_ids" && anyDuplicated(ids$taxID)) {
    # Remove repeated identical taxID/node_id pairs before testing ambiguity.
    mappings <- unique(ids[c("taxID", "node_id")])
    # A source TaxID associated with two different nodes is not a well-defined mapping.
    if (anyDuplicated(mappings$taxID))
      stop("An input taxID maps to multiple native nodes; resolve the ambiguous mapping first.")
  }
  # Reduce the filtered input to one representative row per native node; this defines the
  # unique-node capacity model.
  observations <- ids[!duplicated(ids$node_id), c("taxID", "node_id"), drop = FALSE]
  # Required IDs represent clade coverage constraints, not mandatory selection of those exact
  # internal nodes.
  # Requirements specify clade coverage, not individual sequence records.
  # A required ancestor is redundant when a required descendant exists.
  # Check whether ancestor lies on the path from node toward the root.
  is_ancestor <- function(ancestor, node) {
    # Start at the candidate descendant and move upward through parent links.
    current <- node
    # Keep a per-path visited set to guard against cycles.
    visited <- character()
    # Continue only while the current node is a valid, nonempty ID in the parent lookup.
    while (!is.na(current) && nzchar(current) && current %in% names(parent)) {
      # Reaching ancestor proves that ancestor contains node in its subtree.
      if (identical(current, ancestor)) return(TRUE)
      # A repeated ancestor-walk node indicates a cycle.
      if (current %in% visited) stop("Cycle detected in taxonomy.")
      visited <- c(visited, current)
      # Look up the next ancestor, removing the name attribute.
      next_node <- unname(parent[[current]])
      # A missing parent or a self-parented root terminates the walk.
      if (is.na(next_node) || identical(next_node, current)) break
      current <- next_node
    }
    FALSE
  }
  # Remove ancestor requirements made redundant by a more specific required descendant.
  if (length(required) > 1L) {
    # For each required node, ask whether any OTHER required node is nested beneath it.
    redundant <- vapply(required, function(x) {
      any(vapply(setdiff(required, x), function(y) is_ancestor(x, y),
                 logical(1)))
    }, logical(1))
    # Retain only the nonredundant required clade roots.
    required <- required[!redundant]
  }
  # Compute direct eligible-node indicators; subtree totals are computed below on demand.
  # Postorder counts and required coverage are computed on the observed tree.
  # Count the eligible distinct native-node IDs mapped directly to each observed backbone node.
  direct <- tabulate(match(observations$node_id, nodes$id), nbins = nrow(nodes))
  # Attach native node IDs so direct[node_id] is an accessible named count.
  names(direct) <- nodes$id
  # Read background species counts separately from eligible sampling capacities.
  bg <- stats::setNames(as.numeric(taxlist$spp_df$species_count),
                        as.character(taxlist$spp_df$taxID))
  # Missing background counts are invalid; zero counts are handled later by fallback sampling
  # weights.
  if (anyNA(bg)) stop("Background counts contain NA.")
  # After redundant ancestors are removed, required clades should not nest within each other.
  # Required clades are disjoint after redundant ancestor removal.
  # Count required clade roots inside a given subtree; used to reserve at least one sample per
  # required child clade.
  required_below <- function(node) sum(required %in% descendants(node))
  # Compute capacity as the number of distinct eligible native nodes in a subtree (not
  # background species_count).
  total <- function(node) sum(direct[descendants(node)])
  # Capacity of the selected root determines the maximum without-replacement sample size.
  available <- total(taxon)
  # Reject infeasible requests rather than repeating TaxIDs or silently returning fewer taxa.
  if (!replacement && m > available)
    stop("Requested m (", m, ") exceeds the number of unique eligible taxonomic nodes (",
         available, "). Reduce m or use replacement = TRUE.")
  # Guard against an empty effective sampling request.
  if (m < 1L) stop("No observations available after filtering.")
  # Each disjoint required clade needs at least one selected taxon.
  if (m < length(required)) stop("m is smaller than the number of required clades.")
  # Each required clade must contain at least one eligible native node.
  for (r in required) if (total(r) < 1) stop("Required clade has no observations: ", r)
  # Recursive tree descent chooses native node IDs; selecting source sequences happens only
  # afterward.
  # Recursion selects observed backbone nodes; terminal mapping is separate.
  # draw_node(node, n) returns exactly n sampled native IDs from the subtree rooted at node.
  draw_node <- function(node, n) {
    # A zero budget yields no selected taxa and stops this branch immediately.
    if (n == 0L) return(character())
    # Find children of this node that are inside the permitted root scope.
    kids <- intersect(children[[node]], scope)
    # Drop child subtrees with zero eligible native-node capacity.
    kids <- kids[vapply(kids, total, numeric(1)) > 0]
    # Initially, each eligible child subtree is one allocation option.
    opts <- kids
    # If the internal node itself has an eligible input mapping, add it as a separate direct-
    # sampling option.
    if (direct[node] > 0) opts <- c(opts, node)
    # Fail if the requested branch contains no eligible direct node or descendant.
    if (!length(opts)) stop("No observed taxon below node: ", node)
    # Compute each option’s capacity: direct node count for the current node, subtree count for
    # a child.
    cap <- vapply(opts, function(x) if (x == node) direct[node] else total(x), numeric(1))
    # With replacement, branch allocations are not bounded by the count of unique eligible IDs.
    if (replacement) cap[] <- Inf
    # A required descendant clade forces at least one selection somewhere inside the
    # corresponding child subtree.
    # A required clade is satisfied by any observation in its subtree.
    # Never force a direct observation at an internal required node.
    # The direct current-node option receives no mandatory allocation: a required internal clade
    # can be satisfied by descendants.
    lower <- vapply(opts, function(x) if (x == node) 0L
                    else as.integer(required_below(x)), integer(1))
    # A requirement on the current node is already satisfied if any selected observation comes
    # from this subtree.
    # Required ancestor clades are satisfied by any selected descendant.
    # This branch intentionally does not force the current node itself to be sampled.
    if (node %in% required && !any(lower > 0L)) {
      # A requirement at this node needs one sample anywhere below it.
      # n >= 1 is enforced by the parent allocation.
    }
    # A branch cannot satisfy more mandatory child allocations than its available budget.
    if (n < sum(lower)) stop("Infeasible required allocation below ", node)
    # Select weights from background diversity for known_species, otherwise use uniform weights.
    weights <- if (sampling == "known_species") {
      # A direct observation at this node receives unit weight; child options receive the
      # background species count.
      vapply(opts, function(x) if (x == node) 1 else {
        z <- bg[x]
        # Missing child background weights fall back to zero; nonnegative weights are enforced.
        if (is.na(z)) 0 else max(0, z)
      }, numeric(1))
      # Agnostic mode assigns the same weight to every option.
    } else rep(1, length(opts))
    # Name every vector by its native node option so allocation results can be addressed by node
    # ID.
    names(cap) <- names(lower) <- names(weights) <- opts
    # Allocate the budget among direct-node and child-subtree options, respecting mandatory
    # minima and capacities.
    alloc <- .ts_allocate(n, lower, cap, weights, method, randomize, sampling)
    # Accumulate sampled native IDs from all options.
    result <- character()
    # Only recurse into options assigned a positive number of selections.
    for (x in opts[alloc > 0]) {
      # The direct-node option emits that node ID; repeated copies are possible only when
      # replacement permits them.
      if (x == node) result <- c(result, rep(node, alloc[x]))
      # For a child option, recursively distribute that child’s allocation through its subtree.
      else result <- c(result, draw_node(x, alloc[x]))
    }
    result
  }
  # Start recursion at the user-selected root with the requested budget.
  selected_nodes <- draw_node(taxon, as.integer(m))
  # Check that recursion returned exactly m native IDs.
  if (length(selected_nodes) != m) stop("Internal error: sample size mismatch.")
  # Without replacement, a native node may appear at most once in the output.
  if (!replacement && anyDuplicated(selected_nodes))
    stop("Internal error: duplicate native taxonomic node sampled without replacement.")
  # Only after selecting native nodes do we choose corresponding source TaxIDs and sequence
  # records.
  # Resolve each sampled node to a representative input taxon/sequence.
  # Track which unique native-node representatives remain available for mapping.
  remaining <- observations
  # Preallocate one output record per selected native node.
  output <- vector("list", length(selected_nodes))
  # Map each sampled native node back to a source TaxID and sequence ID.
  for (i in seq_along(selected_nodes)) {
    node <- selected_nodes[i]
    # Find the remaining representative record for this native node.
    choices <- which(remaining$node_id == node)
    # Missing matches indicate an exhausted or inconsistent mapping.
    if (!length(choices)) stop("Internal error: exhausted node ", node)
    # Select a representative index. Because observations was deduplicated by node_id, there
    # should be one choice per node.
    j <- choices[sample.int(length(choices), 1L)]
    # A native node may have several source TaxIDs and sequences; choose one of each only after
    # node selection.
    # Choose one input TaxID for the selected native node, then one sequence.
    # Enumerate unique source TaxIDs associated with this selected native node.
    candidates <- unique(ids$taxID[ids$node_id == node])
    # Choose one source TaxID without interpreting numeric-looking IDs as ranges.
    original <- .ts_draw_one(candidates)
    # Collect all sequence IDs associated with that source TaxID and native node.
    sequence_choices <- ids$seqID[ids$taxID == original & ids$node_id == node]
    # A selected source TaxID must have at least one corresponding sequence.
    if (!length(sequence_choices)) stop("No sequence ID for ", original)
    # Choose a sequence ID from the available mappings.
    seqid <- .ts_draw_one(sequence_choices)
    # Record the selection index and all three identifier levels for traceability.
    output[[i]] <- data.frame(sample_index = i, node_id = node,
                              taxID = original, seqID = seqid,
                              stringsAsFactors = FALSE)
    # Without replacement, remove the selected native-node representative from the remaining
    # mapping pool.
    if (!replacement) remaining <- remaining[-j, , drop = FALSE]
  }
  # Combine all per-selection records into one data.frame.
  output <- do.call(rbind, output)
  # Reset row names to a clean sequential index.
  rownames(output) <- NULL
  # Expose native node IDs as outputNodes and original source TaxIDs as outputIDs.
  taxlist$outputNodes <- output$node_id
  taxlist$outputIDs <- output$taxID
  # Keep the full mapping so callers can inspect the selected node/TaxID/sequence relationships.
  taxlist$outputMapping <- output
  # Store all effective run settings with the result for reproducibility and inspection.
  taxlist$ts.params <- list(taxon = taxon, m = m, method = method,
                            randomize = randomize, replacement = replacement,
                            ignoreIDs = ignored, requireIDs = required,
                            sampling = sampling)
  # Optionally extract selected sequences from the supplied FASTA file.
  if (!is.null(seq_file)) {
    # Sequence extraction uses the selected sequence IDs recorded in outputMapping.
    taxlist <- extract_sequences(taxlist, seq_file, verbose = verbose)
    # If an output FASTA path is supplied, write the extracted sequences using seqinr.
    if (!is.null(out_file)) seqinr::write.fasta(sequences = taxlist$outputSeqs,
                                                names = names(taxlist$outputSeqs),
                                                file.out = out_file)
    # Writing a FASTA file requires an input FASTA file from which sequences can be extracted.
  } else if (!is.null(out_file)) stop("out_file requires seq_file.")
  # Optionally report the number of selected native taxonomic nodes.
  if (verbose) message("Selected ", length(taxlist$outputIDs), " taxonomic nodes.")
  # Return the original taxlist enriched with selections, mapping, parameters, and optional
  # sequences.
  taxlist
}
