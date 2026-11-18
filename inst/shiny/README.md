# TaxonSampling Shiny application (updated)

A standalone Shiny interface for the TaxonSampling R package's DuckDB-based reference backbones and sampling workflow.

## Installation

Install the **current local development version** of TaxonSampling (the supplied `R/` sources, plus its package metadata and dependencies), then install the app dependencies:

```r
install.packages(c("shiny", "shinyFiles", "DT", "ggplot2", "DBI", "duckdb"))
shiny::runApp("TaxonSampling_Shiny_updated")
```

## Workflow

1. Build a reference database from extracted NCBI/OTT taxonomy or a synthetic/custom Newick tree, **or** inspect an existing `.duckdb` file.
2. Upload a headerless, tab-delimited file containing **two columns**: taxon ID and sequence ID. Each taxon can have multiple sequences.
3. Prepare counts using `get_counts()`. The app proposes an observed root node, which you can change.
4. Run `run_TS()` with diversity or experimental balanced sampling and download the sampled pairs or complete R object.
5. Compare input and sampled representation against background diversity, grouped by rank.

Use the **Browse file** and **Browse directory** buttons to choose reference inputs and output locations graphically, or type/paste a path manually. These selectors use `shinyFiles` and browse the filesystem **on the computer running the Shiny server**. For a locally launched application this is your own computer; on a remote server these are server-side paths. Uploaded input files are handled by Shiny. Building a full NCBI/OTT database can require substantial resources and time; prebuilding it outside Shiny is recommended for production.

## Compatibility and caveats

- Based on the source functions uploaded November 13, 2026 (`build_taxonomy_db`, `get_counts`, `run_TS`).
- A two-column input is mandatory in the uploaded `get_counts()` implementation. A one-column ID list is not supported without generated sequence IDs.
- Sampling roots and ignore/require IDs must be **native observed backbone node IDs**, not necessarily the user's original input IDs.
- Balanced sampling can return fewer than `m` and does not support `randomize="after_first_round"`.
- For synthetic trees without OTT species-rank annotations, select `background_unit="tips"`.
- The existing package's `get_summary_df()` and `plot.taxonsampling()` still assume numeric IDs or legacy output semantics; the app implements a backbone-compatible comparison instead.
- Comparisons at a single rank normalize within the *observed* groups at that rank; they do not quantify unsampled clades. `All observed nodes` mixes nested clades, so interpret its frequency comparison cautiously.
- The figure uses log10 relative frequencies and a diagonal. Zero sampled frequencies are omitted from the sampled overlay but remain visible in the table.
- Database creation and sampling are synchronous; large datasets can temporarily block the Shiny session.
- For multi-user deployment, restrict reference build and file-system access to trusted users.

## Original application

The original `ui.R`/`server.R` code is not needed to run this version. The original TS logo is retained.

## Packaging for CRAN

Include `shinyFiles` in the package `DESCRIPTION` `Imports:` alongside `shiny`, `DT`, and `ggplot2`. Place the app under `inst/shiny/TaxonSampling/` and expose it through a package launcher using `system.file()`. The graphical selectors work when launched locally. For remote deployment, filesystem browsing should be restricted or replaced by uploads.

## Interactive pairwise representation plots

The Representation tab provides three independent interactive scatterplots:

1. Full reference backbone vs full input database.
2. Full reference backbone vs sampled dataset.
3. Full input database vs sampled dataset.

Each point represents a taxonomic group at the chosen rank. Hover to inspect taxon
name, identifier, rank, counts, frequencies, and log2 fold change. Orange indicates
overrepresentation on the vertical axis, blue underrepresentation, and grey an
absolute log2 ratio at or below the user-selected cutoff. This is an effect-size
threshold, **not a statistical significance test**. The default threshold of 1
means a twofold difference. A small offset is added only for plotting zeros on
log axes; exact zero frequencies are retained for classification.

Frequencies are calculated from all qualifying groups at the selected rank
before limiting the number of plotted points. `All observed nodes` includes
ancestor-descendant overlap and is not a mutually exclusive partition.

Additional app dependency: `plotly` (as well as `shiny`, `shinyFiles`, `DT`,
`ggplot2`, `DBI`, `duckdb`). For package distribution, add the required
packages to DESCRIPTION Imports.

### Parameter help

Each configurable input has a small **?** icon next to its label. Hover over it (or focus it with the keyboard) to see a short description. These tooltips use native browser titles and require no additional packages.


## Graphical help and count preparation troubleshooting

Hover over or keyboard-focus a `?` icon to see its description. Tooltips are rendered
using CSS and do not require JavaScript or extra R packages.

Before selecting **Prepare sequence counts**, select an existing reference database
in the Backbone tab and click **Inspect database**. Upload a nonempty, headerless,
tab-delimited file containing two columns (taxon ID and sequence ID). For example:

```
9606    human_mt
10090   mouse_mt
```

If counting fails, the app now displays a descriptive error and writes the same
message to the R console. If it still fails, copy the complete console message
and the first few lines of the input file for debugging.


## Hierarchical treemaps

The **Treemaps** tab shows three interactive taxonomic treemaps: full reference, full input sequence database, and sampled sequences. Select a parent rank (e.g. phylum) and a nested rank (e.g. class). The reference plot queries the **complete** DuckDB backbone, including taxa absent from the input. Input and sampled plots use the same backbone groups and their propagated observation counts. Areas encode counts; each plot is scaled independently. Hover shows taxon name, count, and percentage of its parent. Residual diversity not assigned to a child-rank node appears as **Other / unclassified**. The parent/child rank must be ordered correctly. Treemap construction can be slow for large backbones at species level.

The application requires `plotly`, `DBI`, and `duckdb` in addition to the existing dependencies.


### Treemap counts
Treemap input and sampled counts are now propagated over the **complete** reference backbone rather than only nodes retained by `get_counts()`. This prevents classes from being incorrectly placed in Other/unclassified when they are absent from the pruned input hierarchy. The application needs write access to its DuckDB database to create a temporary counting table; the database content is not modified.


## Bias analysis tab

Three interactive plots share the taxonomic-rank, maximum-group, and log2 bias cutoff
settings: ranked original/sample bias bars, improvement in absolute bias, and a
diverging heatmap. Hover for taxon ID, counts, frequencies, and log2 ratios.
Export the full bias table as TSV.

Bias is calculated as log2(dataset frequency / reference frequency). The
improvement metric is abs(original bias) - abs(sampled bias). Positive values
indicate reduced absolute bias. Undefined 0/0 ratios and non-finite improvement
values are not shown as improvements. Infinite bias ratios are capped only for
plotting; exports retain the actual values. The heatmap uses a continuous
blue–grey–orange scale; its grey midpoint indicates zero bias, while the
ranked plot applies the user-defined effect-size cutoff.

These plots reuse the existing comparison-data logic. For the 'All observed
nodes' setting, frequencies include overlapping ancestors; use a named rank
for more interpretable comparisons. R/Shiny runtime testing is required.


### Stable treemap colors

Treemap colors are assigned from taxonomic node identifiers using the full
reference hierarchy, not from abundance. Each parent taxon has a stable hue;
its child taxa have stable related hues. The same node keeps its color in the
reference, full sequence, and sampled treemaps. **Other / unclassified** is
always grey. Rectangle area (not color) represents relative abundance.


### Treemap colors

Parent taxa use a fixed, contrasting Plotly qualitative palette ordered by parent taxon ID. Every child inherits the same color as its parent across all three datasets. Other / unclassified is grey. Rectangle area, not color, represents abundance. When more than ten parents are displayed, palette colors repeat.


### OTT treemap correction (total branch values)

This version preserves `branchvalues="total"` and enforces parent >= sum(children)
for each rendered treemap. Child groups are de-duplicated and, if counts are
inconsistent, offending child rectangles are omitted rather than rescaled or
misrepresented; their parent count remains in **Other / unclassified**.
Warnings identify how many groups were omitted. This is a conservative
visualization safeguard, **not** a repair to source taxonomy counts.
Investigate warnings before using individual child areas for scientific analysis.

### OTT treemap disjoint-clade correction

OTT may contain repeated occurrences of a selected taxonomic rank along an ancestor chain. The treemap now excludes child-rank clades nested inside another child-rank clade under the nearest parent-rank ancestor, avoiding overlapping counts. It uses full background totals, adds an explicit residual, and retains Plotly `branchvalues="total"`. Any remaining inconsistency raises an informative error listing affected parent IDs rather than silently dropping groups.
