# ============================================================================
# TaxonSampling: corrected legacy allocation algorithms
#
# This file implements the allocation policies originally provided by the
# ten legacy sample_* functions.
#
# It is shared by:
#   1. The legacy sampling entry points.
#   2. The modern run_TS() recursive sampling engine.
#
# Its responsibility is to decide HOW MANY taxonomic IDs should be sampled
# from each child branch of a given taxonomic node.
#
# It does not select actual TaxIDs, traverse the taxonomy, or extract sequences.
# Those operations are performed elsewhere.
#
# ---------------------------------------------------------------------------
# INPUT PARAMETERS
# ---------------------------------------------------------------------------
#
# m:
#   Number of ADDITIONAL taxa that must be allocated.
#
# m_i:
#   Named vector containing existing or mandatory allocations for each
#   branch. These allocations are preserved as minimum requirements.
#
# childrenCount:
#   Named vector containing the maximum sampling capacity of each branch.
#   For sampling without replacement, this corresponds to the number of
#   distinct eligible taxonomic IDs in the branch.
#
# childrenCountSpp:
#   Named vector containing background diversity counts for each branch.
#   These values determine sampling probabilities when
#   sampling = "known_species".
#
# method:
#   "diversity" or "balanced".
#
# randomize:
#   "no", "yes", or "after_first_round".
#
# sampling:
#   "agnostic" or "known_species".
#
# ---------------------------------------------------------------------------
# OUTPUT
# ---------------------------------------------------------------------------
#
# A named vector giving the FINAL number of taxa allocated to each branch.
#
# The returned allocation must satisfy:
#
#   allocation[i] >= m_i[i]
#   allocation[i] <= childrenCount[i]
#
# and:
#
#   sum(allocation) = sum(m_i) + m
#
# The original algorithms did not always satisfy these constraints.
# The repaired implementation explicitly enforces them.
# ============================================================================


.ts_legacy_fixed_allocate <- function(
    m,
    m_i,
    childrenCount,
    childrenCountSpp,
    method,
    randomize,
    sampling
) {

  # ========================================================================
  # STEP 1: IDENTIFY THE CHILD BRANCHES
  # ========================================================================

  # Retrieve the names of the child branches.
  #
  # These names are native taxonomy node identifiers supplied by the caller.
  # They allow allocations, capacities, and weights to be matched correctly
  # even when the input vectors have different element orders.

  ids <- names(childrenCount)

  # Each child branch must have a valid, unique identifier.
  #
  # Duplicate or missing names would make it impossible to reliably match
  # allocations to the correct branch.

  if (is.null(ids) || anyDuplicated(ids) || anyNA(ids))
    stop("childrenCount must have unique, non-missing names.")


  # ========================================================================
  # STEP 2: VERIFY THAT ALL INPUT VECTORS DESCRIBE THE SAME BRANCHES
  # ========================================================================

  # Compare branch names across:
  #
  #   m_i              : existing/minimum allocations
  #   childrenCount    : maximum branch capacities
  #   childrenCountSpp : background diversity weights
  #
  # Sorting makes the comparison independent of the original vector order.

  if (!identical(sort(names(m_i)), sort(ids)) ||
      !identical(sort(names(childrenCountSpp)), sort(ids)))
    stop("Allocation vectors must have identical branch names.")


  # ========================================================================
  # STEP 3: STANDARDIZE THE ALLOCATION VECTORS
  # ========================================================================

  # Extract existing/minimum allocations in the order defined by ids.
  #
  # These allocations may represent mandatory coverage of required clades
  # or allocations that have already been reserved by the caller.
  #
  # They must never be reduced by the allocation algorithm.

  lower <- as.integer(m_i[ids])
  names(lower) <- ids

  # Extract maximum capacities in the same branch order.
  #
  # Capacity is distinct from background diversity.
  #
  # Example:
  #
  #   Branch A contains 10,000 species in the reference taxonomy,
  #   but only 15 of those species occur in the input dataset.
  #
  #   Its sampling capacity is 15, not 10,000.
  #
  # The reference diversity may influence its sampling probability,
  # but it cannot increase the number of available input taxa.

  cap <- as.numeric(childrenCount[ids])
  names(cap) <- ids

  # Extract background species counts.
  #
  # These counts are used to weight random allocation decisions in
  # known_species mode.
  #
  # They do not directly determine the maximum allocation capacity.

  bg <- as.numeric(childrenCountSpp[ids])
  names(bg) <- ids


  # ========================================================================
  # STEP 4: VALIDATE ALLOCATION CONSTRAINTS
  # ========================================================================

  # Reject invalid allocation requests.
  #
  # Conditions checked:
  #
  #   - m must be a single nonmissing, nonnegative integer.
  #   - Existing allocations cannot contain missing values.
  #   - Existing allocations cannot be negative.
  #   - Capacities cannot be missing.
  #   - Capacities must accommodate the existing allocations.
  #   - Background weights must be finite and nonnegative.
  #
  # These checks prevent the allocator from operating on inputs that
  # cannot represent a valid sampling problem.

  if (length(m) != 1L || is.na(m) || m < 0 || m != floor(m) ||
      anyNA(lower) || any(lower < 0) || anyNA(cap) || any(cap < lower) ||
      anyNA(bg) || any(!is.finite(bg)) || any(bg < 0))
    stop("Invalid allocation budget, minimum, capacity, or weights.")

  # Check that the additional sampling budget can fit within the remaining
  # capacities of all branches.
  #
  # For each branch:
  #
  #   remaining_capacity[i] = cap[i] - lower[i]
  #
  # Therefore:
  #
  #   m <= sum(cap - lower)
  #
  # must hold for a feasible allocation.
  #
  # Without this check, the original algorithms could return fewer taxa
  # than requested or allocate more taxa than a branch contains.

  if (m > sum(cap - lower))
    stop("Requested allocation exceeds available capacity.")


  # ========================================================================
  # STEP 5: INITIALIZE THE OUTPUT ALLOCATION
  # ========================================================================

  # Start with the mandatory allocations.
  #
  # Every subsequent allocation decision adds taxa to this vector.
  #
  # This guarantees that the initial lower bounds are retained.

  out <- lower

  # If no additional taxa are requested, the existing allocations are
  # already the complete solution.

  if (!m) return(out)


  # ========================================================================
  # STEP 6: DEFINE THE SAMPLING PROBABILITY FUNCTION
  # ========================================================================

  # This helper calculates the probability weights for the currently
  # eligible branches.
  #
  # eligible:
  #   Integer indices identifying branches that can receive another taxon.
  #
  # There are two possible weighting schemes:
  #
  # AGNOSTIC:
  #   Every eligible branch has the same probability.
  #
  # KNOWN_SPECIES:
  #   Branch probabilities are proportional to their background diversity.
  #
  # For known_species mode:
  #
  #                   bg[i]
  #   P(branch i) = -----------
  #                 sum(bg[j])
  #
  # where the denominator includes only currently eligible branches.
  #
  # IMPORTANT:
  #   These are weights for the current draw, not permanent probabilities.
  #   Once a branch reaches capacity, it is removed from future draws.

  prob <- function(eligible) {

    # NULL is interpreted by sample.int() as uniform sampling.
    #
    # Thus, agnostic mode gives equal probability to every eligible branch.

    if (sampling == "agnostic") return(NULL)

    # Extract background diversity counts for eligible branches only.
    #
    # Saturated branches are excluded before this function is called.

    w <- bg[eligible]

    # If every eligible branch has zero background weight, use uniform
    # sampling as a fallback.
    #
    # This avoids an undefined probability distribution when the sum of
    # weights is zero.

    if (sum(w) <= 0) return(NULL)

    # Return the unnormalized weights.
    #
    # sample.int() performs the required normalization internally.

    w
  }


  # ========================================================================
  # STEP 7: DEFINE A HELPER FOR RANDOM BRANCH SELECTION
  # ========================================================================

  # Select one or more eligible branch indices.
  #
  # Parameters:
  #
  # eligible:
  #   Integer indices of branches eligible for selection.
  #
  # size:
  #   Number of branch indices to draw.
  #
  # replace:
  #   Whether the same branch index may be drawn multiple times within
  #   this particular call.
  #
  # Note that replace = FALSE here refers to branch-index selection within
  # one draw. It is not the same as run_TS(replacement = FALSE), which
  # controls repeated selection of taxonomic IDs.
  #
  # This helper also avoids a known R behavior:
  #
  #   sample(5, 1)
  #
  # samples from integers 1:5 rather than returning the value 5.
  #
  # sample.int() operates directly on positional indices, avoiding this
  # ambiguity.

  pick <- function(eligible, size = 1L, replace = FALSE) {

    eligible[
      sample.int(
        length(eligible),
        size,
        replace = replace,
        prob = prob(eligible)
      )
    ]
  }


  # ========================================================================
  # STEP 8: DEFINE A HELPER TO ALLOCATE ONE ADDITIONAL TAXON
  # ========================================================================

  # Randomly choose one eligible branch and increase its allocation by one.
  #
  # This operation is the basic building block of the randomized
  # diversity-first algorithms.
  #
  # The <<- operator modifies out in the enclosing allocation function,
  # rather than creating a new local variable inside add_one().

  add_one <- function(eligible) {

    # Select one eligible branch using the configured weighting scheme.

    j <- pick(eligible)

    # Add one taxon to the selected branch.

    out[j] <<- out[j] + 1L
  }


  # ========================================================================
  # STEP 9: SELECT THE ALLOCATION STRATEGY
  # ========================================================================

  # Two major allocation families are implemented:
  #
  #   1. diversity
  #   2. balanced
  #
  # The following branches implement their respective allocation policies.


  if (method == "diversity") {

    # ======================================================================
    # DIVERSITY-FIRST ALLOCATION
    # ======================================================================
    #
    # The objective is to distribute sampling effort across different
    # lineages before repeatedly sampling from the same lineage.
    #
    # Three modes are available:
    #
    #   randomize = "no"
    #       Complete deterministic coverage rounds followed by a
    #       stochastic remainder.
    #
    #   randomize = "yes"
    #       Sequential stochastic allocation from the first draw.
    #
    #   randomize = "after_first_round"
    #       One complete coverage round, when feasible, followed by
    #       sequential stochastic allocation.
    #
    # These modes reproduce the intended behavior of the original
    # diversity-first algorithms, with corrections for sample-size
    # conservation and branch capacities.


    if (randomize == "no") {

      # ====================================================================
      # STRATEGY: DIVERSITY-FIRST, DETERMINISTIC ROUNDS
      # ====================================================================
      #
      # Corresponds to:
      #
      #   sample_a_d_n()
      #   sample_k_d_n()
      #
      # The algorithm repeatedly assigns one additional taxon to every
      # branch that still has available capacity.
      #
      # Example:
      #
      #   Current allocations: 0, 0, 0, 0
      #   Capacities:          3, 2, 4, 1
      #   Additional budget:   7
      #
      # First round:
      #   1, 1, 1, 1
      #   Remaining budget = 3
      #
      # The fourth branch is now saturated.
      #
      # Since only three branches remain eligible, a second full round
      # can be completed:
      #
      #   2, 2, 2, 1
      #   Remaining budget = 0
      #
      # This creates broad representation before allocating additional
      # observations to already represented branches.

      while (m > 0L) {

        # Identify all branches with remaining capacity.
        #
        # A branch is eligible if:
        #
        #   out[i] < cap[i]

        eligible <- which(out < cap)

        # Stop if every branch has reached its capacity.
        #
        # For a feasible request, this should not occur while a positive
        # budget remains.

        if (!length(eligible)) break

        # A complete round requires one taxon for every eligible branch.
        #
        # If the remaining budget is smaller than the number of eligible
        # branches, stop the complete-round phase and allocate the
        # remainder separately.

        if (m < length(eligible)) break

        # Assign one additional taxon to every eligible branch.

        out[eligible] <- out[eligible] + 1L

        # Subtract the number of allocated taxa from the remaining budget.

        m <- m - length(eligible)
      }

      # ====================================================================
      # ALLOCATE THE REMAINDER
      # ====================================================================
      #
      # If the remaining budget is insufficient for another complete
      # coverage round, choose that many distinct eligible branches.
      #
      # This corrects a defect in the original deterministic
      # diversity-first implementation, which could fail to allocate
      # the remaining observations.

      if (m > 0L) {

        # Recompute eligible branches after all completed rounds.

        eligible <- which(out < cap)

        # Select exactly m distinct eligible branches.
        #
        # The default replace = FALSE ensures that no branch receives
        # more than one additional taxon during this incomplete round.
        #
        # Selection is uniform in agnostic mode or weighted by background
        # diversity in known_species mode.

        chosen <- pick(eligible, size = m)

        # Give one additional taxon to each selected branch.

        out[chosen] <- out[chosen] + 1L
      }

    } else {

      # ====================================================================
      # STOCHASTIC DIVERSITY-FIRST STRATEGIES
      # ====================================================================
      #
      # These include:
      #
      #   randomize = "yes"
      #   randomize = "after_first_round"
      #
      # Both ultimately allocate taxa one at a time using add_one().
      #
      # The difference is whether the algorithm first guarantees one
      # allocation to every eligible branch.


      if (randomize == "after_first_round") {

        # ==================================================================
        # STRATEGY: DIVERSITY-FIRST, RANDOMIZE AFTER FIRST ROUND
        # ==================================================================
        #
        # Corresponds to:
        #
        #   sample_a_d_a()
        #   sample_k_d_a()
        #
        # First, attempt to allocate one taxon to every eligible branch.
        #
        # If the budget is insufficient to cover every branch, skip this
        # guaranteed round and use sequential stochastic allocation.
        #
        # Unlike randomize = "no", subsequent complete coverage rounds
        # are not enforced.

        # Identify branches that can still receive an additional taxon.

        eligible <- which(out < cap)

        # A guaranteed first round is possible only when the remaining
        # budget is at least the number of eligible branches.

        if (m >= length(eligible)) {

          # Give every eligible branch one additional taxon.

          out[eligible] <- out[eligible] + 1L

          # Deduct the number of taxa allocated in this first round.

          m <- m - length(eligible)
        }
      }


      # ====================================================================
      # SEQUENTIAL STOCHASTIC ALLOCATION
      # ====================================================================
      #
      # Corresponds to:
      #
      #   sample_a_d_y()
      #   sample_k_d_y()
      #
      # It also handles the remaining budget for:
      #
      #   sample_a_d_a()
      #   sample_k_d_a()
      #
      # Each iteration:
      #
      #   1. Identifies branches with remaining capacity.
      #   2. Selects one branch according to the probability weights.
      #   3. Adds one taxon to that branch.
      #   4. Reduces the remaining budget by one.
      #
      # Unlike deterministic diversity rounds, a branch may receive
      # multiple allocations before another branch receives its first.
      #
      # Consequently, randomize = "yes" does not guarantee broad
      # coverage at a particular taxonomic level.

      while (m > 0L) {

        # Identify currently eligible branches and allocate one taxon
        # to a randomly selected branch.
        #
        # add_one() updates out in the enclosing function.

        add_one(which(out < cap))

        # One additional taxon has now been allocated.

        m <- m - 1L
      }
    }


  } else if (method == "balanced") {

    # ======================================================================
    # BALANCED ALLOCATION
    # ======================================================================
    #
    # Balanced sampling attempts to distribute the allocation budget
    # across branches.
    #
    # Two modes are implemented:
    #
    #   randomize = "no"
    #       Repeated equal-share allocation with capacity constraints.
    #
    #   randomize = "yes"
    #       Multinomial allocation followed by corrections for mandatory
    #       minima, capacities, and the exact sample-size requirement.
    #
    # IMPORTANT:
    #   The repaired balanced algorithms are not mathematically identical
    #   to the original algorithms in every scenario.
    #
    #   Corrections are necessary because the original implementations
    #   could exceed branch capacities or fail to return the requested
    #   sample size.


    if (randomize == "no") {

      # ====================================================================
      # STRATEGY: DETERMINISTIC BALANCED ALLOCATION
      # ====================================================================
      #
      # Corresponds to:
      #
      #   sample_a_b_n()
      #   sample_k_b_n()
      #
      # The algorithm repeatedly divides the remaining budget as evenly
      # as possible among eligible branches.
      #
      # Branch capacities limit how many taxa each branch can receive.
      #
      # Example:
      #
      #   Capacities:          2, 5, 5
      #   Existing allocation: 0, 0, 0
      #   Additional budget:   9
      #
      # Initial equal share:
      #
      #   floor(9 / 3) = 3
      #
      # Capacity-aware allocation:
      #
      #   2, 3, 3
      #
      # One taxon remains.
      #
      # The first branch is saturated, so the remaining taxon is
      # distributed between the other two branches.

      while (m > 0L) {

        # Identify branches with remaining capacity.

        eligible <- which(out < cap)

        # Stop if all branches are saturated.
        #
        # A feasible allocation should exhaust the budget before
        # exhausting all available capacity.

        if (!length(eligible)) break

        # Calculate the number of taxa that can be assigned equally
        # to every currently eligible branch.
        #
        # Integer division is implemented using floor().
        #
        # For example:
        #
        #   m = 10
        #   eligible branches = 3
        #
        #   whole = floor(10 / 3) = 3

        whole <- floor(m / length(eligible))

        if (whole > 0L) {

          # ================================================================
          # ALLOCATE THE INTEGER SHARE
          # ================================================================
          #
          # The proposed increment is whole for each eligible branch.
          #
          # However, a branch may have fewer than whole available
          # sampling slots.
          #
          # Therefore:
          #
          #   increment[i] = min(
          #     whole,
          #     cap[i] - out[i]
          #   )
          #
          # This prevents allocation above capacity.

          inc <- pmin(whole, cap[eligible] - out[eligible])

          # Add the capacity-constrained increments to the output.

          out[eligible] <- out[eligible] + inc

          # Reduce the remaining budget by the number of taxa
          # actually allocated, rather than the proposed allocation.
          #
          # This is essential when one or more branches become saturated
          # before receiving their full equal share.

          m <- m - sum(inc)

        } else {

          # ================================================================
          # ALLOCATE THE REMAINDER
          # ================================================================
          #
          # This branch executes when:
          #
          #   m < number of eligible branches
          #
          # In that situation, every eligible branch cannot receive
          # another taxon.
          #
          # Select m distinct branches and give each one taxon.

          chosen <- pick(eligible, size = m)

          # Add one taxon to each selected branch.

          out[chosen] <- out[chosen] + 1L

          # The complete remaining budget has been allocated.

          m <- 0L
        }
      }


    } else {

      # ====================================================================
      # STRATEGY: RANDOMIZED BALANCED ALLOCATION
      # ====================================================================
      #
      # Corresponds to:
      #
      #   sample_a_b_y()
      #   sample_k_b_y()
      #
      # This strategy begins with a multinomial allocation proposal.
      #
      # A multinomial draw distributes a fixed number of allocations
      # among several branches according to their probability weights.
      #
      # The original approach is preserved as the initial proposal,
      # but the proposal is subsequently corrected to satisfy:
      #
      #   1. Mandatory minimum allocations.
      #   2. Maximum branch capacities.
      #   3. The exact total sample size.
      #
      # These corrections can alter the original multinomial
      # probability distribution.
      #
      # The direct allocation audit showed that randomized balanced
      # sampling is not statistically equivalent to the original
      # implementation in every tested scenario.


      # ====================================================================
      # STEP A: CALCULATE THE FINAL TARGET ALLOCATION
      # ====================================================================

      # The final sample size includes both:
      #
      #   - Existing mandatory allocations.
      #   - The additional budget m.
      #
      # For example:
      #
      #   lower = c(1, 0, 2)
      #   m = 7
      #
      #   target = 1 + 0 + 2 + 7 = 10

      target <- sum(out) + m


      # ====================================================================
      # STEP B: DEFINE MULTINOMIAL WEIGHTS
      # ====================================================================

      # In known_species mode, use background diversity counts.
      #
      # In agnostic mode, give every branch equal weight.
      #
      # Unlike prob(), this initial multinomial proposal includes
      # all branches, even those whose capacities are already exhausted.
      # Capacity constraints are applied afterward.

      weights <- if (sampling == "known_species") bg else rep(1, length(ids))

      # If every branch has zero background weight, use uniform weights.
      #
      # Otherwise, rmultinom() would not have a valid probability vector.

      if (sum(weights) <= 0) weights[] <- 1


      # ====================================================================
      # STEP C: GENERATE THE MULTINOMIAL PROPOSAL
      # ====================================================================

      # Draw one multinomial allocation vector.
      #
      # The total number of allocations is target.
      #
      # Each branch receives a nonnegative integer count.
      #
      # Before applying constraints:
      #
      #   sum(proposal) = target
      #
      # and, for each branch:
      #
      #   E[proposal[i]] =
      #     target * weights[i] / sum(weights)
      #
      # The multinomial proposal does not automatically respect branch
      # capacities or mandatory allocations.

      proposal <- as.integer(stats::rmultinom(1L, target, weights)[, 1L])

      # Associate each proposed allocation with its branch ID.

      names(proposal) <- ids


      # ====================================================================
      # STEP D: ENFORCE MINIMUM AND MAXIMUM ALLOCATIONS
      # ====================================================================

      # First, limit each proposal to its branch capacity:
      #
      #   min(cap[i], proposal[i])
      #
      # Then enforce the mandatory minimum:
      #
      #   max(lower[i], ...)
      #
      # Since input validation guarantees lower[i] <= cap[i],
      # the resulting allocation satisfies:
      #
      #   lower[i] <= proposal[i] <= cap[i]
      #
      # However, these corrections may change the total allocation.
      #
      # Therefore, additional reconciliation is necessary.

      proposal <- pmax(lower, pmin(cap, proposal))


      # ====================================================================
      # STEP E: REMOVE EXCESS ALLOCATIONS
      # ====================================================================
      #
      # If enforcing mandatory minima increased the proposal total
      # above target, remove allocations until the requested total
      # is restored.
      #
      # Allocations cannot be removed below their mandatory minima.
      #
      # Branches are chosen using the same probability-weighting
      # mechanism defined by pick().
      #
      # NOTE:
      #   This is a stochastic correction procedure.
      #   It does not produce a multinomial distribution conditioned
      #   on the allocation constraints.

      while (sum(proposal) > target) {

        # Only branches above their mandatory minima can lose a taxon.

        eligible <- which(proposal > lower)

        # Select one eligible branch.
        #
        # In agnostic mode, selection is uniform.
        # In known_species mode, selection uses background weights.

        j <- pick(eligible)

        # Remove one allocation from the selected branch.

        proposal[j] <- proposal[j] - 1L
      }


      # ====================================================================
      # STEP F: FILL ANY REMAINING ALLOCATION DEFICIT
      # ====================================================================
      #
      # Capacity clipping can reduce the total below target.
      #
      # In that case, additional taxa are assigned one at a time
      # until the requested sample size is reached.
      #
      # Only branches with remaining capacity are eligible.
      #
      # This preserves capacity constraints while ensuring the
      # final allocation total is correct.

      while (sum(proposal) < target) {

        # Identify branches that can still receive another taxon.

        eligible <- which(proposal < cap)

        # Select one eligible branch using the configured weights.

        j <- pick(eligible)

        # Add one allocation to the selected branch.

        proposal[j] <- proposal[j] + 1L
      }


      # ====================================================================
      # STEP G: STORE THE CORRECTED MULTINOMIAL ALLOCATION
      # ====================================================================

      # Replace the initial minimum-allocation vector with the corrected
      # multinomial allocation.
      #
      # At this point:
      #
      #   sum(proposal) = target
      #
      #   proposal[i] >= lower[i]
      #
      #   proposal[i] <= cap[i]

      out <- proposal
    }

  } else {

    # Reject allocation methods other than "diversity" and "balanced".
    #
    # The public run_TS() interface normally validates this option before
    # calling the allocator, but this check also protects direct calls.

    stop("Unknown allocation method.")
  }


  # ========================================================================
  # STEP 10: RETURN THE FINAL ALLOCATION VECTOR
  # ========================================================================
  #
  # The result contains the number of taxa assigned to each child branch.
  #
  # The recursive sampling engine uses these values to determine how many
  # taxonomic IDs to select from each corresponding subtree.
  #
  # This function does not choose the actual TaxIDs.
  #
  # For feasible inputs and supported allocation modes, the intended
  # postconditions are:
  #
  #   sum(out) = sum(m_i) + original_m
  #
  #   m_i[i] <= out[i] <= childrenCount[i]
  #
  # These constraints are enforced by the allocation procedures above.

  out
}
