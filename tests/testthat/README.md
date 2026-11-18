# TaxonSampling sampling-engine tests

Run from package root after installing dependencies:

```r
install.packages("testthat")
testthat::test_local()
```

Or run `R CMD check .` once the package and dependencies are installed.

The suite uses **in-memory taxonomy fixtures**. No downloaded NCBI/OTT/Synthetic
DuckDB files are needed for the default tests. The fixture tests arbitrary
native IDs, mapping to original TaxIDs, rank-independent clade traversal,
exact sample size, capacities, exclusions, required clades, reproducibility,
allocation statistics, and invalid input.

**Known regression tests**: `test-replacement-regression.R` includes a test
expected to expose the existing direct-ancestor requirement bug. The replacement
test verifies a separate behavior and may pass. These are intentional failures, not assertions
that the current implementation is correct. Do not delete them to get green CI.

Optional performance check (not run by default):

```sh
TS_RUN_PERFORMANCE=1 TS_MAX_SECONDS=30 Rscript -e 'testthat::test_local()'
```

The performance test is a smoke test, not a full scaling benchmark. For a
publication, benchmark on several taxonomy sizes and sample fractions,
recording R version, CPU, memory, elapsed time, and Git commit.

The original `tests/test_cross_identity.R` and
`tests/test_cross_backbone_sampling.R` remain as external integration scripts;
run those separately when full DuckDB reference databases are available.

**Not yet tested by this suite**: end-to-end `get_counts()` mapping with actual
DuckDB fixtures; FASTA extraction; statistical representativeness relative to
an empirical taxonomic background. Those are separate integration and paper
validation tasks.
