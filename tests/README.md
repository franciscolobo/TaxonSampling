# Cross-identity database validation

Run from the TaxonSampling project directory:

```bash
Rscript test_cross_identity.R taxonomy_data/taxonomy_db
```

To also exercise the package's `map_taxa()` API (recommended):

```bash
TS_CROSS_IDENTITY_R=R/build_cross_identity_db.R \
  Rscript test_cross_identity.R taxonomy_data/taxonomy_db
```

If the R files are in a different directory, change `TS_CROSS_IDENTITY_R` accordingly.

Optional second positional argument sets the results directory.

## Tests

- Required schema, nonmissing identifiers, duplicate links.
- Correct stored direct-link degrees and unique/ambiguous classifications.
- NCBI–OTT link provenance checked against `ott.crossreferences` and native node tables.
- OTT–Synthetic links checked against `synthetic.nodes.ott_id`.
- Transitive NCBI–OTT–Synthetic paths checked against the underlying direct links.
- Reports unique species-to-species links separately from all-rank mappings.
- Diagnostic for transitive paths with nonunique endpoints.
- Optional strict `map_taxa()` NCBI→OTT→NCBI round-trip and ambiguity rejection.

Outputs include `test_results.csv`, `direct_mapping_summary.csv`, `transitive_mapping_summary.csv`, species-level coverage CSVs, and `transitive_endpoint_ambiguity.csv`.

**Interpretation:** Tests validate identifier consistency, not equivalence of taxonomic concepts. Mapping coverage is reported without an arbitrary pass/fail threshold. A missing cross-reference does not establish biological absence.

**Performance:** Some tests scan millions of links. Allow time and disk space for DuckDB temporary processing.
