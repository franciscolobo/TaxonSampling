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
