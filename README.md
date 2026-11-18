# TaxonSampling

**Phylogeny-aware taxon sampling for reducing taxonomic bias in biological datasets**

TaxonSampling is an R package designed to select phylogenetically representative subsets of taxa from potentially biased biological datasets.

Public genomic and biological databases contain substantial taxonomic sampling biases. Certain organisms, particularly model species and their close relatives, are disproportionately represented, while large portions of biological diversity remain poorly sampled.

TaxonSampling addresses this problem by using a reference taxonomy or phylogenetic tree to guide the selection of taxa. Rather than sampling organisms uniformly from the available dataset, it considers their distribution across the reference hierarchy, allowing the selected subset to better represent the diversity of the underlying taxonomic or phylogenetic background.

The package supports large-scale taxonomic databases, user-provided phylogenetic trees, alternative sampling strategies, and repeated sampling procedures.

## Key features

- **Phylogeny-aware sampling:** Select taxa while accounting for their distribution across a reference taxonomy or phylogeny.
- **Taxonomic bias reduction:** Reduce the overrepresentation of densely sampled taxonomic groups.
- **Multiple reference backbones:** Support NCBI Taxonomy, Open Tree Taxonomy (OTT), the Open Tree synthetic phylogeny, and custom phylogenetic trees.
- **Scalable database architecture:** Store preprocessed reference hierarchies in persistent DuckDB databases.
- **Reusable reference databases:** Build the taxonomic background once and reuse it across multiple datasets.
- **Flexible sampling strategies:** Support diversity-oriented and balanced sampling approaches.
- **Alternative randomization schemes:** Control the degree of randomization during taxon selection.
- **Sampling with or without replacement:** Enable subsampling and bootstrap-like procedures.
- **Phylogeny-aware resampling:** Generate repeated taxon subsets that better reflect the diversity of the reference background.
- **Sampling diagnostics:** Summarize and visualize the taxonomic representation of sampled datasets.

## Motivation

### Taxonomic bias in biological databases

The availability of genomic and other biological data is strongly influenced by research priorities, sequencing costs, organism accessibility, and historical research practices.

Consequently, the distribution of species represented in biological databases differs substantially from the distribution of known biodiversity.

For example, a dataset containing thousands of species may include hundreds of closely related mammals but comparatively few representatives of many other animal lineages.

Uniform random sampling from this dataset would tend to preserve its original taxonomic biases.

TaxonSampling provides an alternative: use information about the underlying taxonomic or phylogenetic diversity to guide the selection process.

### General principle

Consider a dataset containing 1,000 taxa, distributed across three hypothetical evolutionary lineages:

| Lineage | Available taxa | Reference diversity |
|---|---:|---:|
| A | 800 | 20% |
| B | 150 | 40% |
| C | 50 | 40% |

A uniform random sample of 100 taxa would be expected to contain approximately 80 taxa from lineage A, 15 from lineage B, and 5 from lineage C.

This reproduces the original sampling bias.

A diversity-aware strategy instead uses the reference hierarchy to increase the representation of undersampled lineages, subject to the availability of eligible taxa.

TaxonSampling performs this process hierarchically, using the relationships among taxa rather than relying exclusively on predefined taxonomic ranks.

Importantly, the objective is not necessarily to reproduce the exact proportions of the reference diversity. The objective is to generate samples whose phylogenetic representation is less distorted by the distribution of available observations.

## Supported taxonomic and phylogenetic backbones

TaxonSampling supports four types of reference hierarchies.

| Backbone | Description | Typical applications |
|---|---|---|
| NCBI Taxonomy | Hierarchical taxonomy maintained by NCBI | Genomic datasets indexed by NCBI TaxIDs |
| Open Tree Taxonomy (OTT) | Integrated taxonomy combining multiple taxonomic sources | Broad biodiversity datasets |
| Open Tree synthetic phylogeny | Synthetic tree combining taxonomic and phylogenetic information | Phylogeny-aware sampling across the tree of life |
| Custom tree | User-provided phylogenetic tree | Specialized taxonomic groups, strain-level analyses, and alternative evolutionary hypotheses |

### NCBI Taxonomy

NCBI Taxonomy is particularly useful for genomic applications because most public sequence databases associate their records with NCBI Taxonomy identifiers.

TaxonSampling uses the hierarchical relationships among taxonomic nodes to estimate the distribution of available taxa relative to the reference background.

### Open Tree Taxonomy

The Open Tree Taxonomy integrates information from multiple taxonomic sources into a unified hierarchical classification.

It provides a broader taxonomic framework that is not restricted to organisms represented in NCBI sequence databases.

This can be useful when the objective is to represent biodiversity more generally.

### Open Tree synthetic phylogeny

The Open Tree synthetic phylogeny combines information from published phylogenetic studies and taxonomic relationships.

Unlike a purely taxonomic hierarchy, the synthetic tree incorporates phylogenetic relationships supported by source trees.

TaxonSampling preserves both taxon-associated nodes and synthetic internal nodes when constructing its database representation.

Synthetic internal nodes without taxonomic identifiers are assigned stable internal identifiers.

### Custom phylogenetic trees

Custom trees allow users to apply TaxonSampling beyond existing taxonomic databases.

Potential applications include:

- Sampling strains within microbial species.
- Selecting representatives from gene-family phylogenies.
- Subsampling viral lineages.
- Sampling populations or isolates from an evolutionary tree.
- Evaluating alternative phylogenetic hypotheses.

Custom trees are expected in Newick format.

**Note:** Taxonomic hierarchies and phylogenetic trees are not interchangeable. Taxonomy-based sampling uses classification relationships, whereas tree-based sampling uses the topology of the supplied phylogeny. Neither necessarily provides evolutionary branch lengths or divergence times.

## Installation

TaxonSampling is under active development.

Install the development version from its GitHub repository:

```r
# Install remotes if necessary
install.packages("remotes")

# Install TaxonSampling
remotes::install_github(
  "franciscolobo/TaxonSampling"
)
```

Load the package:

```r
library(TaxonSampling)
```

For large reference databases, sufficient disk space and memory are required during database construction.

The database backend uses DuckDB.

## Overview of the workflow

TaxonSampling separates reference database construction from dataset-specific taxon counting and sampling.

```text
REFERENCE DATA
    |
    +-- NCBI Taxonomy
    |
    +-- Open Tree Taxonomy
    |
    +-- Open Tree synthetic phylogeny
    |
    +-- Custom Newick tree
    |
    v
build_taxonomy_db()
    |
    v
Persistent DuckDB database
    |
    |  Reusable across datasets
    |
    v
get_counts()
    |
    |  User-provided taxon IDs
    |  Mapping and validation
    |  Dataset-specific counts
    |
    v
run_TS()
    |
    |  Sampling strategy
    |  Target sample size
    |  Randomization
    |  Replacement policy
    |
    v
Sampled taxon subset
    |
    +-- Summary statistics
    +-- Diagnostic plots
    +-- Repeated sampling
    +-- Downstream analyses
```

This architecture avoids repeatedly parsing millions of taxonomic nodes when analyzing different datasets.

## 1. Building a reference database

The `build_taxonomy_db()` function processes a reference taxonomy or phylogenetic tree and stores its hierarchical structure in a DuckDB database.

### NCBI Taxonomy

Download and extract the NCBI taxonomy dump from:

https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/

After extracting the taxonomy files:

```r
ncbi_db <- build_taxonomy_db(
  taxonomy_path = "data/ncbi",
  taxonomy = "ncbi",
  output_dir = "taxonomy_db",
  overwrite = FALSE
)
```

This creates:

```text
taxonomy_db/
    ncbi.duckdb
```

### Open Tree Taxonomy

Download and extract an OTT release from the Open Tree of Life project.

```r
ott_db <- build_taxonomy_db(
  taxonomy_path = "data/ott",
  taxonomy = "ott",
  output_dir = "taxonomy_db",
  include_suppressed = TRUE
)
```

The resulting database is stored as:

```text
taxonomy_db/
    ott.duckdb
```

The `include_suppressed` parameter controls the inclusion of suppressed taxa when processing OTT data.

### Open Tree synthetic phylogeny

The synthetic phylogeny can be supplied as a Newick file.

```r
synthetic_db <- build_taxonomy_db(
  taxonomy_path = "data/labelled_supertree_ottnames.tre",
  taxonomy = "synthetic",
  ott_path = "data/ott",
  output_dir = "taxonomy_db"
)
```

The `ott_path` argument supplies OTT metadata for taxonomic annotations, including species-rank identification.

Without OTT metadata, species-rank background counts may be unavailable, although the topology and descendant tip counts can still be represented.

The resulting database is:

```text
taxonomy_db/
    synthetic.duckdb
```

### Custom phylogenetic tree

```r
custom_db <- build_taxonomy_db(
  taxonomy_path = "data/my_tree.nwk",
  taxonomy = "custom",
  output_dir = "taxonomy_db",
  db_name = "my_tree.duckdb"
)
```

Custom trees allow the package to operate on identifiers that are not associated with NCBI or OTT.

## 2. Database structure

The reference databases store the hierarchical relationships and annotations required for sampling.

Principal node attributes include:

| Field | Description |
|---|---|
| `id` | Unique node identifier |
| `parent` | Identifier of the parent node |
| `name` | Taxon or node name |
| `rank` | Taxonomic rank, when available |
| `ott_id` | Associated OTT identifier, when available |
| `node_type` | Classification of the node |
| `label_status` | Information about node-label interpretation |
| `species_count` | Number of descendant species-level taxa |
| `tip_count` | Number of descendant terminal nodes |
| `depth` | Depth of the node within the hierarchy |
| `child_count` | Number of immediate children |

The databases retain both terminal and internal nodes.

In synthetic phylogenies, internal nodes without taxonomic identifiers are preserved rather than discarded.

### Background diversity

Two background measures are particularly important:

**Species counts (`species_count`)**

Represent the number of species-level taxa descending from a node.

This measure is appropriate when species diversity is the relevant reference.

**Tip counts (`tip_count`)**

Represent the number of terminal nodes descending from a node.

This measure can be used for trees containing strains, populations, operational taxonomic units, or other entities without conventional species-rank annotations.

These measures are not necessarily equivalent.

## 3. Preparing the input dataset

TaxonSampling requires a set of taxa eligible for selection.

For example:

```r
taxa <- data.frame(
  taxID = c(
    "9606",
    "9598",
    "10090",
    "10116",
    "7955",
    "7227",
    "6239",
    "3702"
  )
)
```

These identifiers correspond to NCBI Taxonomy IDs.

The input dataset may contain a taxonomically biased collection of organisms, such as species with available genome assemblies, transcriptomes, or protein sequences.

TaxonSampling does not require the input dataset to contain every taxon represented in the reference hierarchy.

## 4. Computing taxonomic counts

The `get_counts()` function associates input taxa with a previously constructed reference database and computes the information required for sampling.

```r
count_data <- get_counts(
  taxonomy_db = "taxonomy_db/ncbi.duckdb",
  ids_df = taxa,
  id_type = "ncbi",
  count_mode = "input_ids",
  background_unit = "species"
)
```

The function separates the observed distribution of input taxa from the background diversity represented by the reference hierarchy.

### Identifier types

The `id_type` argument specifies the identifier namespace:

| Value | Description |
|---|---|
| `"ncbi"` | NCBI Taxonomy identifiers |
| `"ott"` | Open Tree Taxonomy identifiers |
| `"native"` | Identifiers native to the selected backbone |

### Counting modes

The `count_mode` argument controls the units used for the input-side counting procedure.

- `"input_ids"`: Counts based on the supplied input identifiers.
- `"backbone_nodes"`: Counts based on the corresponding reference backbone nodes.

### Background units

The `background_unit` argument selects the reference diversity measure:

- `"species"`: Descendant species-level taxa.
- `"tips"`: Descendant terminal nodes.

For conventional species-level genomic datasets, `"species"` is generally the appropriate background measure.

For custom trees without taxonomic ranks, `"tips"` may be more appropriate.

### Returned object

The resulting object contains the information needed by the sampling algorithm.

Its principal components include:

```r
names(count_data)
```

```text
nodes
ids_df
countIDs
spp_df
mapping
unmapped_ids
metadata
```

These components provide the hierarchical structure, dataset-specific counts, background information, and identifier-mapping results.

## 5. Running TaxonSampling

The main sampling function is `run_TS()`.

It uses the information prepared by `get_counts()` to select a subset of the available taxa.

The sampling procedure can be configured using several options.

### Sampling methods

**Diversity-oriented sampling**

```r
method = "diversity"
```

Prioritizes the representation of taxonomic or phylogenetic diversity using the reference hierarchy.

This is the primary strategy for reducing taxonomic sampling bias.

**Balanced sampling**

```r
method = "balanced"
```

Provides an alternative strategy intended to balance representation across hierarchical groups.

The balanced strategy is experimental and may return fewer taxa than requested under some conditions.

For applications requiring a fixed sample size, the diversity-oriented strategy is generally preferable.

### Randomization

The `randomize` argument supports:

| Value | Behavior |
|---|---|
| `"no"` | No randomization |
| `"yes"` | Randomized selection |
| `"after_first_round"` | Randomization after the initial sampling round |

Randomized sampling is particularly useful for generating repeated representative subsets.

### Sampling modes

The `sampling` argument supports:

- `"agnostic"`
- `"known_species"`

These modes determine how species-related information is incorporated into the sampling procedure.

### Replacement

The `replacement` argument controls whether previously selected taxa can be selected again.

Sampling without replacement is useful when constructing a nonredundant representative dataset.

Sampling with replacement can support bootstrap-like analyses.

### Example workflow

```r
library(TaxonSampling)

# Load eligible taxa
taxa <- data.frame(
  taxID = c(
    "9606", "9598", "10090", "10116",
    "7955", "7227", "6239", "3702"
  )
)

# Prepare reference counts
count_data <- get_counts(
  taxonomy_db = "taxonomy_db/ncbi.duckdb",
  ids_df = taxa,
  id_type = "ncbi",
  count_mode = "input_ids",
  background_unit = "species"
)

# Inspect the prepared data
names(count_data)

# See the installed version's sampling arguments
args(run_TS)
```

The exact `run_TS()` call should be constructed using the sampling parameters and the installed version's function signature.

## 6. Summarizing sampling results

TaxonSampling provides methods for examining the selected taxa and the distribution of the resulting sample.

```r
# Summarize sampling results
summary(result)

# Extract summary data
summary_df <- get_summary_df(result)

# Plot results
plot(result)
```

The resulting summaries can be used to compare the reference diversity, the original dataset, and the selected subset.

### Evaluating taxonomic representation

One useful diagnostic is to compare the relative frequency of taxa belonging to different taxonomic groups before and after sampling.

For a group \(i\), define:

\[
f_i^{\mathrm{input}} =
\frac{n_i^{\mathrm{input}}}{N^{\mathrm{input}}}
\]

\[
f_i^{\mathrm{sample}} =
\frac{n_i^{\mathrm{sample}}}{N^{\mathrm{sample}}}
\]

\[
f_i^{\mathrm{background}} =
\frac{n_i^{\mathrm{background}}}{N^{\mathrm{background}}}
\]

A useful visualization compares the observed frequencies against the reference frequencies.

For instance, scatterplots can display log-transformed frequencies with taxon labels and a diagonal representing equal frequencies.

A representative sample should generally show reduced distortion relative to the reference distribution, although exact correspondence is not guaranteed.

## 7. Applications

### 7.1 Representative genome selection

Genome databases contain many closely related species, especially within intensively studied lineages.

TaxonSampling can help select smaller sets of genomes that retain broader evolutionary representation.

Potential applications include:

- Comparative genomics.
- Gene-family evolution.
- Orthology inference.
- Protein-family characterization.
- Genome annotation benchmarking.
- Training datasets for comparative sequence models.

### 7.2 Phylogeny-aware bootstrap and resampling

Traditional bootstrap procedures sample observations without explicitly accounting for their evolutionary relationships.

When biological observations are phylogenetically related, uniform resampling may repeatedly favor densely represented lineages.

TaxonSampling can be used to generate repeated samples whose taxonomic or phylogenetic representation is informed by a broader reference hierarchy.

A general procedure is:

1. Define the complete set of eligible taxa.
2. Select an appropriate reference taxonomy or phylogeny.
3. Generate repeated diversity-aware samples.
4. Perform the downstream analysis independently for each sample.
5. Evaluate the stability of the results across samples.

Potential applications include assessing the robustness of genotype–phenotype associations, comparative genomic statistics, evolutionary trait analyses, and machine-learning models.

**Statistical consideration:** Phylogeny-aware resampling is not automatically equivalent to conventional nonparametric bootstrapping. Its inferential properties depend on the sampling design, the target population, and the statistic being estimated. Confidence intervals and hypothesis tests require appropriate validation.

### 7.3 Improving phylogenetic marker datasets

Phylogenetic marker collections and genome completeness benchmarks can be affected by the taxonomic composition of their reference genomes.

For example, BUSCO relies on lineage-specific sets of conserved orthologous genes.

TaxonSampling could help construct more representative initial genome collections for developing or evaluating such datasets.

Potential benefits include reducing the dominance of heavily sequenced lineages and evaluating how marker selection changes across different representative genome subsets.

This is a proposed application requiring empirical validation.

### 7.4 Comparative genotype–phenotype analyses

Comparative analyses frequently combine observations from species with shared evolutionary histories.

TaxonSampling can help investigate whether results depend disproportionately on taxonomically overrepresented groups.

Repeated sampling can be used to evaluate the stability of associations across alternative subsets of species.

TaxonSampling complements, but does not replace, explicit phylogenetic comparative methods such as phylogenetic generalized least squares.

### 7.5 Machine-learning dataset construction

Biological machine-learning datasets often contain phylogenetically redundant observations.

This can produce training sets dominated by closely related organisms and can contribute to overly optimistic evaluations when related taxa are split across training and test sets.

TaxonSampling can be used to construct phylogenetically diverse training datasets or to investigate model sensitivity to taxonomic representation.

Additional phylogeny-aware data partitioning is still necessary to prevent evolutionary information leakage between training and evaluation datasets.

## 8. Validation

The effectiveness of TaxonSampling should be evaluated by comparing its outputs against alternative sampling strategies.

Possible baselines include:

- Uniform random sampling.
- Taxonomic stratified sampling.
- Phylogenetic diversity maximization.
- Other representative-selection algorithms.

### Suggested evaluation metrics

**Sample-size accuracy**

Whether the algorithm returns the requested number of taxa.

**Reference distribution similarity**

How closely the sampled taxonomic distribution resembles the selected background distribution.

**Taxonomic coverage**

The number of represented taxonomic groups at different hierarchical levels.

**Phylogenetic coverage**

The distribution of selected taxa across the reference topology.

**Sampling reproducibility**

The variability of results across repeated runs.

**Computational performance**

Runtime, memory consumption, and scalability with reference database size and input dataset size.

No single metric captures all aspects of representative sampling. Evaluation should account for the tradeoff between taxonomic representation, sample size, and available taxa.

## 9. Scalability

TaxonSampling uses persistent DuckDB databases to avoid repeatedly reconstructing large reference hierarchies.

Reference preprocessing and dataset-specific sampling are separate operations.

This design is particularly important when working with millions of taxonomic nodes.

The current database architecture has been exercised with reference hierarchies containing millions of nodes, including NCBI Taxonomy, OTT, and the Open Tree synthetic phylogeny.

Database construction is typically more expensive than subsequent dataset-specific counting and sampling.

Performance depends on the selected backbone, hardware, input dataset, and sampling configuration.

## 10. Limitations

### Reference hierarchy quality

Sampling results depend on the accuracy and completeness of the selected taxonomy or phylogeny.

Different reference backbones may produce different samples.

### Taxonomy versus phylogeny

Taxonomic relationships do not necessarily reflect phylogenetic branch lengths or divergence times.

Consequently, taxonomic representation should not automatically be interpreted as maximization of evolutionary distance.

### Availability constraints

TaxonSampling can only select taxa present in the eligible input dataset.

If a lineage is poorly represented or absent from the input, the algorithm cannot recover its missing diversity.

### Species definitions

Species boundaries and taxonomic ranks differ among organisms and databases.

Background species counts should therefore be interpreted in the context of the selected reference.

### Identifier mapping

Mapping taxa across NCBI, OTT, and synthetic phylogenetic backbones may be incomplete or ambiguous.

Unmapped identifiers should be examined before interpreting sampling results.

### Statistical interpretation

A representative sample does not automatically remove phylogenetic non-independence or all forms of ascertainment bias.

Downstream analyses should continue to use appropriate statistical models.

## 11. Development roadmap

TaxonSampling is undergoing active development.

Planned improvements include:

- Further validation of sampling strategies across reference backbones.
- Improved handling of ambiguous cross-database taxonomic mappings.
- Additional support for custom phylogenetic trees.
- Expanded visualization and diagnostic tools.
- Repeated sampling and phylogeny-aware bootstrap workflows.
- Benchmarking against alternative taxon-selection algorithms.
- Improved computational performance for very large datasets.
- Expanded automated testing and reproducible examples.

## 12. Data sources and references

### NCBI Taxonomy

Schoch CL et al. (2020). NCBI Taxonomy: a comprehensive update on curation, resources and tools. *Database*, baaa062.

https://doi.org/10.1093/database/baaa062

NCBI Taxonomy download:

https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/

### Open Tree of Life

Hinchliff CE et al. (2015). Synthesis of phylogeny and taxonomy into a comprehensive tree of life. *Proceedings of the National Academy of Sciences*, 112, 12764–12769.

https://doi.org/10.1073/pnas.1423041112

Rees JA & Cranston K (2017). Automated assembly of a reference taxonomy for phylogenetic data synthesis. *Biodiversity Data Journal*, 5, e12581.

https://doi.org/10.3897/BDJ.5.e12581

Open Tree of Life:

https://tree.opentreeoflife.org/

Reference taxonomy repository:

https://github.com/OpenTreeOfLife/reference-taxonomy

### DuckDB

https://duckdb.org/

## 13. Citation

If you use TaxonSampling in your research, please cite the software repository and the associated publication when available.

A formal software citation will be provided upon publication.

## 14. Contributing

Contributions, bug reports, and suggestions are welcome through the project's GitHub repository.

Particularly valuable contributions include:

- Reproducible examples of taxonomic sampling biases.
- Benchmark datasets.
- Alternative sampling strategies.
- Performance improvements.
- Support for additional taxonomic databases.
- Validation of phylogeny-aware resampling procedures.

## 15. License

See the `LICENSE` file in the repository for the applicable licensing terms.

---

**TaxonSampling: selecting fewer taxa while representing more of life's diversity.**
