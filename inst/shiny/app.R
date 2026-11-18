library(shiny)
library(TaxonSampling)
library(DT)
library(ggplot2)
library(plotly)
library(shinyFiles)
source("diversity_metrics.R", local = TRUE)

# Accessible parameter labels with compact native-browser hover help.
# No extra JavaScript or Shiny dependencies are required.
parameter_label <- function(label, description) {
  tags$span(
    label, " ",
    tags$span(
      "?", `data-help` = description, `aria-label` = paste(label, description),
      tabindex = "0", role = "note", class = "ts-help-icon"
    )
  )
}

# Persist only local paths, never database contents, in the user's R config directory.
.ts_config_file <- function() {
  dir <- tools::R_user_dir("TaxonSampling", which = "config")
  file.path(dir, "shiny_paths.rds")
}
.ts_read_paths <- function() {
  path <- .ts_config_file()
  if (!file.exists(path)) return(list())
  tryCatch({
    x <- readRDS(path)
    if (!is.list(x)) return(list())
    x
  }, error = function(e) list())
}
.ts_save_paths <- function(workspace = NULL, database = NULL) {
  path <- .ts_config_file()
  x <- .ts_read_paths()
  if (!is.null(workspace) && nzchar(trimws(workspace))) {
    x$workspace <- normalizePath(path.expand(workspace), mustWork = FALSE)
  }
  if (!is.null(database) && nzchar(trimws(database)) && file.exists(database)) {
    x$database <- normalizePath(database, mustWork = TRUE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = "shiny_paths_", tmpdir = dirname(path))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("Could not save Shiny database paths: ", path)
  invisible(x)
}
.ts_saved_paths <- .ts_read_paths()
.ts_default_workspace <- if (is.character(.ts_saved_paths$workspace) &&
    length(.ts_saved_paths$workspace) == 1L &&
    dir.exists(.ts_saved_paths$workspace)) .ts_saved_paths$workspace else "taxonomy_data"
.ts_default_database <- if (is.character(.ts_saved_paths$database) &&
    length(.ts_saved_paths$database) == 1L &&
    file.exists(.ts_saved_paths$database)) .ts_saved_paths$database else
    file.path(.ts_default_workspace, "taxonomy_db", "ncbi.duckdb")

# Paths to local reference files are entered on the host running Shiny.
ui <- fluidPage(
  titlePanel("TaxonSampling | phylogeny-aware sampling"),
      tags$head(tags$style(HTML("
        .ts-help-icon { display:inline-flex; align-items:center; justify-content:center;
          width:15px; height:15px; margin-left:3px; border:1px solid #607d91;
          border-radius:50%; color:#37627c; font-size:11px; font-weight:700;
          line-height:1; cursor:help; vertical-align:middle; background:#f5f8fa; }
        .ts-help-icon { position:relative; }
        .ts-help-icon:hover, .ts-help-icon:focus { color:#fff; background:#37627c; outline-offset:2px; }
        .ts-help-icon:hover::after, .ts-help-icon:focus::after {
          content:attr(data-help); position:absolute; z-index:99999;
          top:calc(100% + 7px); left:0; width:max-content; max-width:300px;
          min-width:180px; white-space:normal; box-sizing:border-box;
          background:#263747; color:#fff; padding:8px 10px;
          border-radius:5px; box-shadow:0 3px 12px rgba(0,0,0,.25);
          font-size:12px; line-height:1.45; font-weight:400;
          text-align:left; pointer-events:none;
        }
      ")),
      tags$script(HTML("Shiny.addCustomMessageHandler('ts_setup_progress',function(x){var a=document.getElementById('ts-activity-live');var l=document.getElementById('ts-log-live');if(a){a.textContent=x.text;if(x.busy!==undefined){var sp=document.getElementById('ts-spinner-live');if(sp)sp.style.display=x.busy?'inline-block':'none';}}if(l&&x.append){var line=document.createElement('div');line.textContent=x.append;l.appendChild(line);l.scrollTop=l.scrollHeight;}});")),
      tags$style(HTML(".ts-working{display:flex;align-items:center;gap:9px;margin:8px 0;font-weight:600}.ts-spinner{display:inline-block;width:19px;height:19px;border:3px solid #d5e1e9;border-top-color:#2376a5;border-radius:50%;animation:ts-spin .7s linear infinite}@keyframes ts-spin{to{transform:rotate(360deg)}}"))),
  sidebarLayout(
    sidebarPanel(width = 4,
      tags$img(src = "TS_logo.png", style = "max-width:180px;max-height:100px;"),
      tabsetPanel(id = "workflow",
        tabPanel("1. Taxonomy databases",
          p("Download raw NCBI, OTT, and OpenTree Synthetic sources and build the three taxonomy databases plus their cross-identity database."),
          textInput("setup_dir", parameter_label("Taxonomy workspace (server path)", "Raw files are stored in raw/ and databases in taxonomy_db/ under this directory."),  .ts_default_workspace),
          shinyDirButton("browse_setup_dir", "Browse workspace", "Choose taxonomy workspace"),
          checkboxGroupInput("setup_sources", "Taxonomy sources", choices=c("NCBI Taxonomy"="ncbi", "Open Tree Taxonomy (OTT)"="ott", "OpenTree Synthetic phylogeny"="synthetic"), selected=c("ncbi", "ott", "synthetic")),
          checkboxInput("setup_download", "Download raw files before building", TRUE),
          checkboxInput("setup_cross", "Build cross-identity database (fourth database)", TRUE),
          checkboxInput("include_suppressed", parameter_label("Include suppressed OTT taxa", "Include flagged OTT taxa in background species counts."), TRUE),
          checkboxInput("overwrite", parameter_label("Overwrite existing files/databases", "Re-download and rebuild selected resources; potentially expensive."), FALSE),
          actionButton("setup_all", "Download and build databases", class="btn-primary"),
          actionButton("refresh_db_stats", "Refresh database statistics"),
          tags$hr(),
          selectInput("active_backbone", "Sampling backbone", choices=c("NCBI"="ncbi", "OTT"="ott", "Synthetic"="synthetic")),
          actionButton("load_backbone", "Load selected backbone"),
          tags$hr(),
          tags$details(
            tags$summary("Advanced: custom tree or existing database"),
            textInput("taxonomy_path", "Extracted taxonomy directory / Newick file", ""),
            shinyFilesButton("browse_taxonomy_file", "Browse file", "Select Newick or taxonomy file", multiple=FALSE),
            shinyDirButton("browse_taxonomy_dir", "Browse directory", "Select extracted taxonomy directory"),
            selectInput("taxonomy", "Reference type", c("NCBI"="ncbi", "OTT"="ott", "Synthetic"="synthetic", "Custom Newick"="custom")),
            textInput("ott_path", "OTT directory for synthetic annotations (optional)", ""),
            shinyDirButton("browse_ott_dir", "Browse OTT directory", "Select OTT directory"),
            textInput("output_dir", "Advanced database output directory", "taxonomy_db"),
            shinyDirButton("browse_output_dir", "Browse output directory", "Select database output directory"),
            textInput("db_name", "Database filename (optional)", ""),
            actionButton("build", "Build single database"),
            textInput("db_path", "Existing DuckDB file",  .ts_default_database),
            shinyFilesButton("browse_db_file", "Browse database", "Select DuckDB database", multiple=FALSE),
            actionButton("inspect", "Inspect existing database")
          )
        ),
        tabPanel("2. Input taxa",
          fileInput("taxa_file", parameter_label("Taxon / sequence ID TSV (headerless, two columns)", "Upload a headerless file with two columns: taxon identifier and sequence/observation identifier. Multiple rows may share a taxon."), accept=c(".tsv", ".txt", ".csv")),
          helpText("One row per taxon–sequence pair. IDs are strings; retain leading zeros. Use the same taxon with multiple sequence IDs if needed."),
          selectInput("id_type", parameter_label("Identifier namespace", "Identifier system used by the uploaded taxa: NCBI TaxIDs, OTT IDs, or identifiers native to the selected backbone."), c("NCBI"="ncbi", "OTT"="ott", "Native backbone ID"="native")),
          selectInput("count_mode", parameter_label("Observation unit", "Choose whether input counts represent distinct supplied taxon IDs or distinct mapped nodes in the reference backbone."), c("Distinct input taxon IDs"="input_ids", "Distinct mapped backbone nodes"="backbone_nodes")),
          selectInput("background_unit", parameter_label("Background unit", "Reference diversity measure: descendant species-level taxa or terminal tips. Tips are useful for custom trees."), c("Species"="species", "Terminal tips"="tips")),
          selectInput("missing_ids", parameter_label("Unmapped identifiers", "What to do when an uploaded identifier cannot be mapped to the reference: warn, drop, or stop with an error."), c("Warn"="warn", "Drop"="drop", "Error"="error")),
          selectInput("duplicate_ids", parameter_label("Repeated identical taxon/sequence pairs", "How to handle duplicate rows with the same taxon and sequence identifiers."), c("Keep first"="first", "Error"="error")),
          actionButton("count", "Prepare dataset counts"),
          downloadButton("download_counts", "Download count object (.rds)")
        ),
        tabPanel("3. Sampling",
          textInput("root", parameter_label("Root node ID (native backbone identifier)", "Root of the clade to sample from, using the native node ID in the observed nodes table."), ""),
          helpText("Use an ID from the observed nodes table. OTT IDs in synthetic trees may use the ott: prefix."),
          numericInput("m", parameter_label("Number of observations", "Requested number of taxon/sequence observations to select. Some strategies may return fewer."), 100, min=1, step=1),
          selectInput("method", parameter_label("Method", "Diversity prioritizes representation across the reference hierarchy; balanced is experimental and may undersample."), c("Diversity"="diversity", "Balanced (experimental)"="balanced")),
          selectInput("randomize", parameter_label("Randomization", "No: deterministic selection; Yes: randomized selection; After first round: randomize subsequent selection rounds."), c("After first round"="after_first_round", "Yes"="yes", "No"="no")),
          selectInput("sampling", parameter_label("Sampling mode", "Known species uses species-level information; agnostic uses the hierarchy without requiring species-level interpretation."), c("Known species"="known_species", "Agnostic"="agnostic")),
          checkboxInput("replacement", parameter_label("Sample with replacement", "Allow the same eligible observation to be selected more than once."), FALSE),
          textAreaInput("ignore", parameter_label("Exclude clade IDs (comma / newline separated)", "Native IDs of clades to exclude from sampling; separate multiple IDs with commas or new lines."), "", rows=2),
          textAreaInput("require", parameter_label("Require clade IDs (comma / newline separated)", "Native IDs of clades to require during sampling; separate multiple IDs with commas or new lines."), "", rows=2),
          numericInput("seed", parameter_label("Random seed", "Seed used to make randomized sampling reproducible when inputs and software versions are unchanged."), 1, min=0, step=1),
          actionButton("sample", "Run TaxonSampling"),
          downloadButton("download_selected", "Download sampled taxa (.tsv)"),
          downloadButton("download_result", "Download full result (.rds)")
        )

      )
    ),
    mainPanel(width=8,
      verbatimTextOutput("status"),
      tabsetPanel(
        tabPanel("Activity log",
          p("Taxonomy setup, input-taxon mapping, and sampling operations."),
          tags$div(class="ts-working", role="status", `aria-live`="polite",
            tags$span(id="ts-spinner-live", class="ts-spinner", style="display:none", `aria-hidden`="true"),
            tags$span(id="ts-activity-live", "Idle")),
          tags$div(style="max-height:480px;overflow-y:auto;white-space:pre-wrap;word-break:break-word;background:#f7f9fa;border:1px solid #ddd;padding:10px;font-family:monospace;font-size:12px;",
            tags$div(id="ts-log-live"))
        ),
        tabPanel("Database statistics",
          tabsetPanel(id = "database_statistics_tabs",
            tabPanel("Taxonomy databases",
              p("Statistics describe each complete taxonomy database, not the subset of input taxa."),
              DTOutput("database_statistics"),
              h4("Cross-taxonomy identity mappings"),
              p("Mapping counts and ambiguity statistics; not phylogenetic node statistics."),
              DTOutput("cross_identity_statistics")),
            tabPanel("Reference taxonomy",
              h4("Reference taxonomy by rank"),
              p("Number of nodes and terminal tips at each taxonomic rank in the selected reference database."),
              DTOutput("reference_table"))
          )
        ),
        tabPanel("Input taxa",
          tabsetPanel(id = "input_taxa_tabs",
            tabPanel("Observed nodes",
              p("Nodes in the induced taxonomy and their propagated counts."),
              DTOutput("nodes_table")),
            tabPanel("Mappings",
              p("Mapping of input identifiers to nodes in the selected reference taxonomy."),
              DTOutput("mapping_table"))
          )
        ),
        tabPanel("Sampling analysis",
          tabsetPanel(id = "sampling_analysis_tabs",
            tabPanel("Selected taxa", DTOutput("selected_table")),
            tabPanel("Representation",
          h4("Representation controls"),
          p("These settings also affect Bias analysis and Diversity metrics where applicable."),
          selectInput("level", parameter_label("Taxonomic rank", "Taxonomic rank used to aggregate and compare frequencies; All observed nodes includes overlapping ancestral groups."), c("Phylum"="phylum", "Class"="class", "Order"="order", "Family"="family", "Genus"="genus", "Species"="species", "All observed nodes"="all")),
          numericInput("top_n", parameter_label("Maximum plotted groups (0 = all)", "Maximum number of taxa shown in each comparison, ranked by input abundance. Set to zero to plot all groups."), 40, min=0, max=100000),
          numericInput("bias_cutoff", parameter_label("Bias cutoff: absolute log2 fold change", "Minimum absolute log2 frequency ratio to color a group as over- or underrepresented; 1 means a twofold difference."), 1, min=0, step=0.25),
          helpText("Orange = overrepresented (Y versus X); blue = underrepresented; grey = within cutoff. Cutoff 1 = twofold. This is an effect-size threshold, not a statistical significance test."),
          checkboxInput("labels", parameter_label("Show taxon labels", "Display taxon names next to plotted points; hover over points for full details even when labels are hidden."), FALSE),
          downloadButton("download_summary", "Download comparison table (.tsv)"),
          tags$hr(),
          h4("1. Full reference vs full input database"), plotlyOutput("comparison_reference_input", height="490px"),
          h4("2. Full reference vs sampled dataset"), plotlyOutput("comparison_reference_sample", height="490px"),
          h4("3. Full input database vs sampled dataset"), plotlyOutput("comparison_input_sample", height="490px"),
          DTOutput("summary_table")),
        tabPanel("Bias analysis",
          h4("Original and sampled taxonomic bias"),
          p("Bias is log2(dataset frequency / reference frequency). Positive values indicate overrepresentation. Negative values indicate underrepresentation."),
          p("All three views use the taxonomic rank, group limit, and bias cutoff in Sampling analysis / Representation."),
          h4("1. Ranked bias comparison"),
          plotlyOutput("bias_ranked", height="600px"),
          h4("2. Change in absolute bias after sampling"),
          plotlyOutput("bias_improvement", height="560px"),
          h4("3. Taxonomic bias heatmap"),
          plotlyOutput("bias_heatmap", height="560px"),
          downloadButton("download_bias", "Download bias analysis (.tsv)"),
          tags$hr(),
          p("Bias uses frequency ratios, not a hypothesis test. Groups with zero frequency have infinite log2 bias; a finite display cap is used for plotting, while the downloaded table retains the original values. Comparisons at 'All observed nodes' include overlapping ancestors and are not independent taxonomic partitions.")
        ),
        tabPanel("Diversity metrics",
          h4("Taxonomic diversity and reference representativeness"),
          p("Metrics use the complete set of groups at the rank selected in Sampling analysis / Representation; the maximum plotted groups setting does not affect calculations."),
          p("Hill q=0 measures represented groups; q=1 is the effective Shannon diversity; q=2 emphasizes dominant groups. JSD (base 2; 0–1) measures distance from the reference frequency distribution."),
          DTOutput("diversity_metrics_table"),
          downloadButton("download_diversity_metrics", "Download metrics (.tsv)"),
          h4("Effective diversity (Hill q=1)"),
          plotlyOutput("diversity_hill_plot", height="380px"),
          h4("Divergence from the reference"),
          plotlyOutput("diversity_jsd_plot", height="380px"),
          verbatimTextOutput("diversity_interpretation"),
          tags$hr(),
          p("Interpretation: frequencies are normalized among groups present at the chosen rank, not over the entire taxonomy. Missing/unranked lineages outside that rank are excluded. 'All observed nodes' cannot be used because overlapping ancestors do not form a partition. Species-level or sparsely populated ranks may give misleading comparisons; use a shared rank with adequate coverage. Counts may represent distinct IDs or backbone nodes according to the input settings, rather than sequence records. Reference group coverage is the fraction of reference groups represented, not a sample-size-adjusted richness estimate.")
        ),
        tabPanel("Treemaps",
          h4("Hierarchical taxonomic representation"),
          p("Choose a parent rank and its nested child rank. Rectangle areas represent counts; hover to inspect taxa."),
          fluidRow(
            column(6, selectInput("tree_root_rank", parameter_label("Parent rank", "Top-level taxonomic groups displayed as outer rectangles."), c("Phylum"="phylum", "Class"="class", "Order"="order", "Family"="family", "Genus"="genus"), selected="phylum")),
            column(6, selectInput("tree_child_rank", parameter_label("Nested rank", "Descendant taxonomic rank nested within each parent group. Must be lower than the parent rank."), c("Class"="class", "Order"="order", "Family"="family", "Genus"="genus", "Species"="species"), selected="class"))
          ),
          p("The full reference uses the complete DuckDB backbone; input and sample use their mapped counts. Only parent-rank lineages not covered by the displayed child-rank groups appear as Other / unclassified. Child-rank groups include all their descendants."),
          h4("1. Full reference database"), plotlyOutput("tree_reference", height="520px"),
          h4("2. Full sequence database"), plotlyOutput("tree_input", height="520px"),
          h4("3. Sampled sequences"), plotlyOutput("tree_sampled", height="520px")
        ),
          )
        ),
        tabPanel("Run details", verbatimTextOutput("details"))
      )
    )
  )
)

server <- function(input, output, session) {
  # Graphical filesystem selectors. These browse the host running Shiny.
  # Use local runApp() for desktop use; remote hosting exposes server paths.
  roots <- shinyFiles::getVolumes()()
  home <- path.expand("~")
  if (dir.exists(home)) roots <- c(Home=home, roots)

  shinyFiles::shinyFileChoose(input, "browse_taxonomy_file", roots=roots,
                              session=session, filetypes=c("tre", "tree", "nwk", "newick", "txt", "tsv", "dmp"))
  shinyFiles::shinyDirChoose(input, "browse_taxonomy_dir", roots=roots,
                             session=session)
  shinyFiles::shinyDirChoose(input, "browse_ott_dir", roots=roots,
                             session=session)
  shinyFiles::shinyDirChoose(input, "browse_output_dir", roots=roots,
                             session=session)
  shinyFiles::shinyDirChoose(input, "browse_setup_dir", roots=roots, session=session)
  shinyFiles::shinyFileChoose(input, "browse_db_file", roots=roots,
                              session=session, filetypes=c("duckdb", "db"))

  observeEvent(input$browse_taxonomy_file, {
    selected <- shinyFiles::parseFilePaths(roots, input$browse_taxonomy_file)
    if (nrow(selected)) updateTextInput(session, "taxonomy_path", value=as.character(selected$datapath[1]))
  })
  observeEvent(input$browse_db_file, {
    selected <- shinyFiles::parseFilePaths(roots, input$browse_db_file)
    if (nrow(selected)) updateTextInput(session, "db_path", value=as.character(selected$datapath[1]))
  })
  bind_directory <- function(selector, field) {
    observeEvent(input[[selector]], {
      selected <- shinyFiles::parseDirPath(roots, input[[selector]])
      if (length(selected) && nzchar(selected[1]))
        updateTextInput(session, field, value=as.character(selected[1]))
    }, ignoreNULL=TRUE)
  }
  bind_directory("browse_taxonomy_dir", "taxonomy_path")
  bind_directory("browse_ott_dir", "ott_path")
  bind_directory("browse_output_dir", "output_dir")
  bind_directory("browse_setup_dir", "setup_dir")
  state <- reactiveValues(db=NULL, counts=NULL, result=NULL, message="Choose an existing reference database or build one.")
  # Restore the last selected reference database when it is still available.
  observeEvent(TRUE, {
    saved <- .ts_read_paths()
    db <- saved$database
    if (is.character(db) && length(db) == 1L && file.exists(db)) {
      tryCatch({
        con <- DBI::dbConnect(duckdb::duckdb(), dbdir=db, read_only=TRUE)
        valid <- tryCatch(all(c("nodes", "metadata") %in% DBI::dbListTables(con)),
                          finally=DBI::dbDisconnect(con, shutdown=TRUE))
        if (valid) {
          state$db <- normalizePath(db, mustWork=TRUE)
          state$message <- paste("Restored reference database:", state$db)
        }
      }, error=function(e) {
        state$message <- paste("Saved database could not be opened:", conditionMessage(e))
      })
    }
  }, once=TRUE)
  report_error <- function(e) {
    msg <- conditionMessage(e)
    if (!nzchar(trimws(msg))) {
      msg <- paste0("Operation interrupted (", paste(class(e), collapse=", "), "). ",
                    "Check that a reference database is loaded and an input file is selected.")
    }
    state$message <- paste("Error:", msg)
    message("[TaxonSampling Shiny] ", state$message)
    if (inherits(e, "shiny.silent.error")) {
      message("[TaxonSampling Shiny] A Shiny req()/validation interrupted the operation.")
    }
    showNotification(state$message, type="error", duration=NULL)
  }
  parse_ids <- function(s) { x <- trimws(unlist(strsplit(s, "[,;\\n\\r\\t ]+"))); unique(x[nzchar(x)]) }
  # A setup workspace always uses setup_taxonsampling()'s documented layout:
  #   workspace/raw/{ncbi,ott,opentree}
  #   workspace/taxonomy_db/{ncbi,ott,synthetic,cross_identity}.duckdb
  db_stats_version <- reactiveVal(0L)
  setup_results <- reactiveVal(NULL)
  log_setup <- function(text, busy=TRUE) {
    session$sendCustomMessage("ts_setup_progress", list(
      text=text, busy=busy,
      append=paste0(format(Sys.time(), "%H:%M:%S"), "  ", text)))
  }
  db_stats <- reactive({
    db_stats_version()
    workspace <- trimws(input$setup_dir)
    validate(need(nzchar(workspace), "Specify a taxonomy workspace."))
    db_dir <- file.path(workspace, "taxonomy_db")
    types <- c("ncbi", "ott", "synthetic")
    do.call(rbind, lapply(types, function(type) {
      path <- file.path(db_dir, paste0(type, ".duckdb"))
      row <- data.frame(database=type, status=if (file.exists(path)) "Available" else "Missing",
                        total_nodes=NA_real_, terminal_nodes=NA_real_,
                        species_at_root=NA_real_, tips_at_root=NA_real_,
                        maximum_depth=NA_real_, root_id=NA_character_,
                        build_date=NA_character_, size_mb=NA_real_,
                        path=normalizePath(path, mustWork=FALSE), stringsAsFactors=FALSE)
      if (!file.exists(path)) return(row)
      row$size_mb <- round(as.numeric(file.info(path)$size) / 1024^2, 2)
      tryCatch({
        con <- DBI::dbConnect(duckdb::duckdb(), dbdir=path, read_only=TRUE)
        on.exit(DBI::dbDisconnect(con, shutdown=TRUE), add=TRUE)
        tables <- DBI::dbListTables(con)
        if ("metadata" %in% tables) {
          meta <- DBI::dbGetQuery(con, "SELECT key, value FROM metadata")
          if (all(c("key", "value") %in% names(meta))) {
            get_meta <- function(k) { v <- meta$value[meta$key == k]; if (length(v)) as.character(v[[1L]]) else NA_character_ }
            row$build_date <- get_meta("build_date")
            row$root_id <- get_meta("root_id")
          }
        }
        if ("nodes" %in% tables) {
          vals <- DBI::dbGetQuery(con, paste(
            "SELECT COUNT(*) AS total_nodes,",
            "SUM(CASE WHEN child_count = 0 THEN 1 ELSE 0 END) AS terminal_nodes,",
            "MAX(depth) AS maximum_depth FROM nodes"))
          row$total_nodes <- as.numeric(vals$total_nodes[[1L]])
          row$terminal_nodes <- as.numeric(vals$terminal_nodes[[1L]])
          row$maximum_depth <- as.numeric(vals$maximum_depth[[1L]])
          if (!is.na(row$root_id) && nzchar(row$root_id)) {
            roots <- DBI::dbGetQuery(con,
              "SELECT species_count, tip_count FROM nodes WHERE id = ? LIMIT 1",
              params=list(row$root_id))
            if (nrow(roots)) {
              row$species_at_root <- as.numeric(roots$species_count[[1L]])
              row$tips_at_root <- as.numeric(roots$tip_count[[1L]])
            }
          }
        }
        row
      }, error=function(e) { row$status <- paste("Inspection error:", conditionMessage(e)); row })
    }))
  })
  output$database_statistics <- renderDT({
    DT::datatable(db_stats(), rownames=FALSE,
                  options=list(pageLength=4, scrollX=TRUE, dom="tip"))
  })
  cross_identity_stats <- reactive({
    db_stats_version()
    workspace <- trimws(input$setup_dir)
    validate(need(nzchar(workspace), "Specify a taxonomy workspace."))
    path <- file.path(workspace, "taxonomy_db", "cross_identity.duckdb")
    metric <- c("Database status", "File size (MB)", "Build date",
                "NCBI → OTT links", "OTT → Synthetic links",
                "Unique direct links", "Ambiguous direct links",
                "NCBI → Synthetic paths", "Unique transitive paths",
                "Ambiguous transitive paths")
    value <- rep(NA_character_, length(metric))
    value[1L] <- if (file.exists(path)) "Available" else "Missing"
    if (!file.exists(path)) return(data.frame(metric, value))
    value[2L] <- as.character(round(as.numeric(file.info(path)$size) / 1024^2, 2))
    tryCatch({
      con <- DBI::dbConnect(duckdb::duckdb(), dbdir=path, read_only=TRUE)
      on.exit(DBI::dbDisconnect(con, shutdown=TRUE), add=TRUE)
      tables <- DBI::dbListTables(con)
      if ("metadata" %in% tables) {
        meta <- DBI::dbGetQuery(con, "SELECT key, value FROM metadata")
        date <- meta$value[meta$key == "build_date"]
        if (length(date)) value[3L] <- as.character(date[[1L]])
      }
      if ("identity_links_status" %in% tables) {
        counts <- DBI::dbGetQuery(con, paste(
          "SELECT",
          "COUNT(*) FILTER (WHERE source_a='ncbi' AND source_b='ott') AS ncbi_ott,",
          "COUNT(*) FILTER (WHERE source_a='ott' AND source_b='synthetic') AS ott_synthetic,",
          "COUNT(*) FILTER (WHERE mapping_status='unique') AS unique_links,",
          "COUNT(*) FILTER (WHERE mapping_status='ambiguous') AS ambiguous_links",
          "FROM identity_links_status"))
        value[4:7] <- as.character(unlist(counts[1L, ], use.names=FALSE))
      }
      if ("transitive_links" %in% tables) {
        counts <- DBI::dbGetQuery(con, paste(
          "SELECT COUNT(*) AS paths,",
          "COUNT(*) FILTER (WHERE mapping_status='unique_path') AS unique_paths,",
          "COUNT(*) FILTER (WHERE mapping_status='ambiguous_path') AS ambiguous_paths",
          "FROM transitive_links"))
        value[8:10] <- as.character(unlist(counts[1L, ], use.names=FALSE))
      }
      data.frame(metric, value)
    }, error=function(e) {
      value[1L] <- paste("Inspection error:", conditionMessage(e))
      data.frame(metric, value)
    })
  })
  output$cross_identity_statistics <- renderDT({
    DT::datatable(cross_identity_stats(), rownames=FALSE,
                  colnames=c("Statistic", "Value"),
                  options=list(dom="t", paging=FALSE, ordering=FALSE))
  })
  observeEvent(input$refresh_db_stats, {
    db_stats_version(db_stats_version() + 1L)
  })
  observeEvent(input$setup_all, {
    session$sendCustomMessage("ts_setup_progress", list(text="Starting database setup...", busy=TRUE))
    on.exit(session$sendCustomMessage("ts_setup_progress", list(
      text="Database setup finished (see activity log)", busy=FALSE)), add=TRUE)
    tryCatch({
      workspace <- trimws(input$setup_dir)
      if (!nzchar(workspace)) stop("Specify a taxonomy workspace.")
      sources <- input$setup_sources
      if (!length(sources)) stop("Select at least one taxonomy source.")
      state$message <- "Preparing taxonomy databases."
      shiny::withProgress(message="TaxonSampling database setup", detail="Starting...", value=0, {
        steps <- length(sources) * (if (isTRUE(input$setup_download)) 3L else 1L) +
          length(sources) + as.integer(isTRUE(input$setup_cross)) + 1L
        completed <- 0L
        update_progress <- function(text) {
          log_setup(text)
          completed <<- completed + 1L
          shiny::setProgress(value=min(0.98, completed / max(steps, 1L)), detail=text)
        }
        result <- TaxonSampling::setup_taxonsampling(
          path=workspace, sources=sources,
          download=isTRUE(input$setup_download),
          build_cross_identity=isTRUE(input$setup_cross),
          overwrite=isTRUE(input$overwrite),
          include_suppressed=isTRUE(input$include_suppressed), verbose=TRUE,
          progress_callback=update_progress)
        setup_results(result)
        .ts_save_paths(workspace=workspace)
        db_stats_version(db_stats_version() + 1L)
        state$message <- paste("Setup finished:", paste(paste(result$database, result$status, sep="="), collapse="; "))
        log_setup(state$message)
        shiny::setProgress(value=1, detail="Finished")
        showNotification(state$message, type=if (any(result$status %in% c("failed", "unavailable", "missing_source", "missing_tree"))) "warning" else "message", duration=12)
      })
    }, error=function(e) { log_setup(paste("ERROR:", conditionMessage(e))); report_error(e) })
  })
  observeEvent(input$load_backbone, {
    tryCatch({
      path <- file.path(trimws(input$setup_dir), "taxonomy_db", paste0(input$active_backbone, ".duckdb"))
      if (!file.exists(path)) stop("Database not found: ", path, ". Build it first.")
      con <- DBI::dbConnect(duckdb::duckdb(), dbdir=path, read_only=TRUE)
      valid <- all(c("nodes", "metadata") %in% DBI::dbListTables(con))
      DBI::dbDisconnect(con, shutdown=TRUE)
      if (!valid) stop("Selected file is not a TaxonSampling taxonomy database.")
      state$db <- normalizePath(path, mustWork=TRUE)
      .ts_save_paths(workspace=trimws(input$setup_dir), database=state$db)
      state$counts <- NULL; state$result <- NULL
      updateTextInput(session, "db_path", value=state$db)
      state$message <- paste("Sampling backbone loaded:", state$db)
    }, error=report_error)
  })
  observeEvent(input$build, {
    tryCatch({
      req(nzchar(trimws(input$taxonomy_path)))
      path <- TaxonSampling::build_taxonomy_db(
        taxonomy_path=trimws(input$taxonomy_path), taxonomy=input$taxonomy,
        output_dir=trimws(input$output_dir),
        db_name=if(nzchar(trimws(input$db_name))) trimws(input$db_name) else NULL,
        ott_path=if(input$taxonomy == "synthetic" && nzchar(trimws(input$ott_path))) trimws(input$ott_path) else NULL,
        include_suppressed=isTRUE(input$include_suppressed), overwrite=isTRUE(input$overwrite)
      )
      state$db <- as.character(path)
      .ts_save_paths(database=state$db)
      updateTextInput(session, "db_path", value=state$db)
      state$counts <- NULL; state$result <- NULL
      state$message <- paste("Reference database created:", state$db)
    }, error=report_error)
  })
  observeEvent(input$inspect, {
    tryCatch({
      path <- normalizePath(trimws(input$db_path), mustWork=TRUE)
      con <- DBI::dbConnect(duckdb::duckdb(), dbdir=path, read_only=FALSE)
      on.exit(DBI::dbDisconnect(con, shutdown=TRUE))
      if (!all(c("nodes", "metadata") %in% DBI::dbListTables(con))) stop("Not a TaxonSampling reference database")
      state$db <- path; state$counts <- NULL; state$result <- NULL
      .ts_save_paths(database=state$db)
      state$message <- paste("Reference database loaded:", path)
    }, error=report_error)
  })
  reference <- reactive({
    req(state$db)
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir=state$db, read_only=FALSE)
    on.exit(DBI::dbDisconnect(con, shutdown=TRUE))
    DBI::dbGetQuery(con, "SELECT rank, COUNT(*) AS nodes, SUM(CASE WHEN node_type = 'tip' THEN 1 ELSE 0 END) AS tips FROM nodes GROUP BY rank ORDER BY nodes DESC LIMIT 100")
  })
  output$reference_table <- renderDT(datatable(reference(), options=list(pageLength=15, scrollX=TRUE)))
  observeEvent(input$count, {
    session$sendCustomMessage("ts_setup_progress", list(
      text="Loading input taxa...", busy=TRUE,
      append=paste0(format(Sys.time(), "%H:%M:%S"), "  Loading input taxa...")))
    on.exit(session$sendCustomMessage("ts_setup_progress", list(
      text="Input-taxon processing finished", busy=FALSE)), add=TRUE)
    tryCatch({
      if (is.null(state$db) || !length(state$db) || !file.exists(state$db)) {
        stop("No reference database loaded. In Taxonomy databases, load a sampling backbone first.", call.=FALSE)
      }
      if (is.null(input$taxa_file) || !nrow(input$taxa_file)) {
        stop("No sequence database uploaded. Select a two-column headerless TSV file first.", call.=FALSE)
      }
      uploaded_path <- input$taxa_file$datapath[[1]]
      if (!file.exists(uploaded_path)) {
        stop("Uploaded file is unavailable. Please upload it again.", call.=FALSE)
      }
      # Fail early with a useful message if an older installed package is loaded.
      required_args <- c("taxonomy_db", "ids_file", "id_type", "count_mode",
                         "background_unit", "missing_ids", "duplicate_ids")
      installed_args <- names(formals(TaxonSampling::get_counts))
      missing_args <- setdiff(required_args, installed_args)
      if (length(missing_args)) {
        stop(paste0(
          "The installed TaxonSampling::get_counts() is incompatible with this Shiny app. ",
          "Missing parameters: ", paste(missing_args, collapse=", "), ". ",
          "Reinstall the current local TaxonSampling package, restart R, and verify ",
          "args(TaxonSampling::get_counts). Installed package path: ",
          find.package("TaxonSampling")
        ), call.=FALSE)
      }
      # File upload path is transient; get_counts reads it immediately.
      # Forward backend messages and warnings to the right-side Activity log.
      counts <- withCallingHandlers(
        TaxonSampling::get_counts(
          taxonomy_db=state$db, ids_file=uploaded_path,
          id_type=input$id_type, count_mode=input$count_mode,
          background_unit=input$background_unit, missing_ids=input$missing_ids,
          duplicate_ids=input$duplicate_ids, verbose=TRUE
        ),
        message=function(cnd) {
          msg <- trimws(conditionMessage(cnd))
          if (nzchar(msg)) log_setup(msg)
          invokeRestart("muffleMessage")
        },
        warning=function(cnd) {
          msg <- trimws(conditionMessage(cnd))
          if (nzchar(msg)) log_setup(paste0("Warning: ", msg))
          invokeRestart("muffleWarning")
        }
      )
      if (!is.list(counts) || is.null(counts$nodes) || is.null(counts$ids_df)) {
        stop("get_counts() returned an unexpected result; verify the installed TaxonSampling version.", call.=FALSE)
      }
      state$counts <- counts; state$result <- NULL
      roots <- counts$nodes$id[is.na(counts$nodes$parent) | !(counts$nodes$parent %in% counts$nodes$id) | counts$nodes$parent == counts$nodes$id]
      if (length(roots)) updateTextInput(session, "root", value=as.character(roots[1]))
      state$message <- sprintf("Counts prepared: %s mapped taxon/sequence pairs; %s observed nodes; %s unmapped IDs.", nrow(counts$ids_df), nrow(counts$nodes), length(counts$unmapped_ids))
      log_setup(state$message)
    }, error=function(e) {
      log_setup(paste("ERROR:", conditionMessage(e)))
      report_error(e)
    })
  })
  output$nodes_table <- renderDT({req(state$counts); datatable(state$counts$nodes, options=list(pageLength=12, scrollX=TRUE), filter="top")})
  output$mapping_table <- renderDT({req(state$counts); datatable(state$counts$mapping, options=list(pageLength=12, scrollX=TRUE))})
  observeEvent(input$sample, {
    tryCatch({
      req(state$counts)
      if (input$method == "balanced" && input$randomize == "after_first_round") stop("Balanced sampling does not support 'after_first_round'.")
      root <- trimws(input$root)
      if (!nzchar(root)) stop("Enter a root node ID")
      if (!root %in% state$counts$nodes$id) stop("Root must be an ID in the observed nodes table")
      set.seed(as.integer(input$seed))
      result <- TaxonSampling::run_TS(
        taxlist=state$counts, taxon=root, m=as.integer(input$m),
        method=input$method, randomize=input$randomize,
        replacement=isTRUE(input$replacement), ignoreIDs=parse_ids(input$ignore),
        requireIDs=parse_ids(input$require), sampling=input$sampling
      )
      state$result <- result
      state$message <- sprintf("Sampling complete: selected %d observations (requested %d).", nrow(result$outputMapping), as.integer(input$m))
    }, error=report_error)
  })
  output$selected_table <- renderDT({req(state$result); datatable(state$result$outputMapping, options=list(pageLength=15, scrollX=TRUE))})
  output$status <- renderText(state$message)
  output$details <- renderPrint({
    if (!is.null(state$counts)) {
      cat("Reference:", state$counts$metadata$taxonomy, "\n")
      cat("Count mode:", state$counts$metadata$count_mode, "\n")
      cat("Unmapped IDs:", paste(head(state$counts$unmapped_ids, 30), collapse=", "), "\n")
    }
    if (!is.null(state$result)) print(state$result$ts.params)
  })
  # Counts are already propagated to ancestors by get_counts(); sampling
  # counts must instead be propagated from selected backbone node IDs.
  comparison_data <- reactive({
    req(state$counts)
    x <- state$counts
    nodes <- x$nodes
    nodes$id <- as.character(nodes$id); nodes$parent <- as.character(nodes$parent)
    bg <- x$spp_df$species_count[match(nodes$id, as.character(x$spp_df$taxID))]
    obs <- as.numeric(x$countIDs[match(nodes$id, names(x$countIDs))]); obs[is.na(obs)] <- 0
    sampled <- numeric(nrow(nodes))
    if (!is.null(state$result)) {
      sampled <- as.numeric(table(factor(state$result$outputNodes, levels=nodes$id)))
      parent_index <- match(nodes$parent, nodes$id)
      depth <- nodes$depth
      for (i in order(depth, decreasing=TRUE)) {
        p <- parent_index[i]
        if (!is.na(p) && p != i) sampled[p] <- sampled[p] + sampled[i]
      }
    }
    data.frame(id=nodes$id, name=nodes$name, rank=nodes$rank,
               background=as.numeric(bg), input=obs, sampled=sampled,
               stringsAsFactors=FALSE)
  })
  filtered_comparison <- reactive({
    d <- comparison_data()
    if (input$level != "all") {
      d <- d[!is.na(d$rank) & tolower(d$rank) == input$level, , drop=FALSE]
    }
    d <- d[is.finite(d$background) & d$background >= 0 &
             is.finite(d$input) & d$input >= 0 &
             is.finite(d$sampled) & d$sampled >= 0, , drop=FALSE]
    if (!nrow(d)) return(d)
    # At a named rank, compare proportions within that rank. For "all",
    # overlapping ancestor/descendant groups are NOT independent partitions.
    totals <- colSums(d[, c("background", "input", "sampled"), drop=FALSE])
    for (nm in c("background", "input", "sampled")) {
      d[[paste0(nm, "_freq")]] <- if (totals[[nm]] > 0) d[[nm]] / totals[[nm]] else 0
    }
    d[order(d$input, decreasing=TRUE), , drop=FALSE]
  })
  output$summary_table <- renderDT({
    d <- filtered_comparison()
    datatable(d, options=list(pageLength=12, scrollX=TRUE), filter="top")
  })

  diversity_metrics_data <- reactive({
    req(state$result)
    validate(need(input$level != "all",
      "Choose a specific taxonomic rank in Sampling analysis / Representation. Overlapping ancestral nodes cannot be used for diversity metrics."))
    d <- filtered_comparison()
    validate(need(nrow(d) > 0, "No groups at this taxonomic rank."))
    validate(need(all(colSums(d[, c("background", "input", "sampled")]) > 0),
      "All three datasets must have positive counts at the selected rank."))
    ts_diversity_metrics(d, input$level)
  })
  output$diversity_metrics_table <- renderDT({
    datatable(diversity_metrics_data(), rownames=FALSE,
      options=list(dom='t', scrollX=TRUE),
      caption="Diversity and divergence, calculated on the full selected-rank distribution")
  })
  output$diversity_hill_plot <- plotly::renderPlotly({
    d <- diversity_metrics_data()
    plotly::plot_ly(d, x=~dataset, y=~hill_q1, type='bar',
      marker=list(color=c('#636EFA','#EF553B','#00CC96')),
      text=~paste0('Dataset: ',dataset,'<br>Hill q=0: ',signif(hill_q0,5),
                   '<br>Hill q=1: ',signif(hill_q1,5),
                   '<br>Hill q=2: ',signif(hill_q2,5)), hoverinfo='text') |>
      plotly::layout(yaxis=list(title='Effective number of taxonomic groups'),
                     xaxis=list(title=''), showlegend=FALSE)
  })
  output$diversity_jsd_plot <- plotly::renderPlotly({
    d <- diversity_metrics_data()
    plotly::plot_ly(d, x=~dataset, y=~jsd_vs_reference, type='bar',
      marker=list(color=c('#636EFA','#EF553B','#00CC96')),
      text=~paste0('JSD versus reference: ',signif(jsd_vs_reference,5)),
      hoverinfo='text') |>
      plotly::layout(yaxis=list(title='Jensen–Shannon divergence (bits)', range=c(0,1)),
                     xaxis=list(title=''), showlegend=FALSE)
  })
  output$diversity_interpretation <- renderText({
    d <- diversity_metrics_data()
    reduction <- d$jsd_reduction_pct[3]
    if (is.na(reduction))
      'JSD reduction is undefined because the original database already matches the reference (JSD = 0).'
    else sprintf('Relative JSD reduction after sampling: %.2f%% (positive indicates closer agreement with the reference).', reduction)
  })
  output$download_diversity_metrics <- downloadHandler(
    filename=function() paste0('taxonsampling_diversity_', input$level, '.tsv'),
    content=function(file) write.table(diversity_metrics_data(), file, sep='\t',
                                      quote=FALSE, row.names=FALSE, na='NA')
  )

  # Each plot uses its own reference (X) and comparison (Y) distribution.
  # The offset is for rendering zeros only; bias classification uses raw
  # frequencies, so zeros are correctly treated as complete absence.
  make_comparison <- function(x_key, y_key, x_title, y_title) {
    d <- filtered_comparison()
    validate(need(nrow(d) > 0, "No taxa at this rank. Select another rank."))
    validate(need(sum(d[[x_key]]) > 0 && sum(d[[y_key]]) > 0,
                  "Both datasets must contain observations at this rank. Run sampling first for sampled comparisons."))
    nmax <- input$top_n
    if (!is.null(nmax) && nmax > 0 && nrow(d) > nmax) {
      # Stable subset across all three plots: based on input abundance.
      d <- head(d, nmax)
    }
    fx <- d[[paste0(x_key, "_freq")]]
    fy <- d[[paste0(y_key, "_freq")]]
    positive <- c(fx[fx > 0], fy[fy > 0])
    offset <- if (length(positive)) min(positive) / 10 else 1e-12
    offset <- max(offset, .Machine$double.xmin)
    log_ratio <- ifelse(fx == 0 & fy == 0, 0,
                 ifelse(fx == 0, Inf, ifelse(fy == 0, -Inf, log2(fy / fx))))
    cutoff <- input$bias_cutoff
    status <- ifelse(log_ratio > cutoff, "Overrepresented",
              ifelse(log_ratio < -cutoff, "Underrepresented", "Within cutoff"))
    status[fx == 0 & fy == 0] <- "Absent in both"
    status[status == "Absent in both"] <- "Within cutoff"
    d$classification <- factor(status, levels=c("Overrepresented", "Underrepresented", "Within cutoff"))
    d$plot_x <- log10(fx + offset)
    d$plot_y <- log10(fy + offset)
    d$tooltip <- paste0(
      "Taxon: ", d$name, "<br>ID: ", d$id,
      "<br>Rank: ", ifelse(is.na(d$rank), "unranked", d$rank),
      "<br>", x_title, " count: ", format(d[[x_key]], big.mark=",", scientific=FALSE),
      "<br>", x_title, " frequency: ", signif(fx, 5),
      "<br>", y_title, " count: ", format(d[[y_key]], big.mark=",", scientific=FALSE),
      "<br>", y_title, " frequency: ", signif(fy, 5),
      "<br>log2(Y/X): ", ifelse(is.infinite(log_ratio), ifelse(log_ratio > 0, "+Inf", "-Inf"), signif(log_ratio, 4)),
      "<br>Classification: ", status)
    # Consistent coordinate limits guarantee a true y=x diagonal.
    bounds <- range(c(d$plot_x, d$plot_y), finite=TRUE)
    if (diff(bounds) == 0) bounds <- bounds + c(-0.5, 0.5)
    p <- ggplot(d, aes(x=plot_x, y=plot_y, text=tooltip, color=classification)) +
      geom_abline(slope=1, intercept=0, linetype="dashed", color="grey45") +
      geom_point(size=2.8, alpha=0.85) +
      scale_color_manual(values=c("Overrepresented"="#E68624", "Underrepresented"="#2479B5", "Within cutoff"="#9B9B9B"), drop=FALSE) +
      coord_equal(xlim=bounds, ylim=bounds) +
      labs(x=paste0("log10(", x_title, " frequency + offset)"),
           y=paste0("log10(", y_title, " frequency + offset)"),
           color="Representation",
           subtitle=paste0("|log2(Y/X)| cutoff: ", cutoff,
                           "   •   zero plotting offset: ", signif(offset, 3))) +
      theme_minimal(base_size=12) + theme(legend.position="bottom")
    if (isTRUE(input$labels)) {
      p <- p + geom_text(aes(label=name), check_overlap=TRUE, hjust=-0.1,
                         size=2.8, show.legend=FALSE)
    }
    plotly::ggplotly(p, tooltip="text")
  }
  output$comparison_reference_input <- plotly::renderPlotly({
    make_comparison("background", "input", "Full reference", "Full input")
  })
  output$comparison_reference_sample <- plotly::renderPlotly({
    req(state$result)
    make_comparison("background", "sampled", "Full reference", "Sampled")
  })
  output$comparison_input_sample <- plotly::renderPlotly({
    req(state$result)
    make_comparison("input", "sampled", "Full input", "Sampled")
  })
  # Shared bias-analysis table: same frequency denominators as scatterplots.
  # Non-finite ratios are retained in the export, and capped ONLY for rendering.
  bias_analysis <- reactive({
    req(state$result)
    d <- filtered_comparison()
    validate(need(nrow(d) > 0, "No taxonomic groups available at this rank."))
    validate(need(sum(d$background_freq) > 0 && sum(d$input_freq) > 0 &&
                  sum(d$sampled_freq) > 0,
                  "Reference, input, and sampled data must all have observations at this rank."))
    ratio <- function(y, x) {
      out <- rep(NA_real_, length(x))
      both <- x == 0 & y == 0
      out[both] <- NA_real_ # Undefined, not zero.
      out[x == 0 & y > 0] <- Inf
      out[y == 0 & x > 0] <- -Inf
      valid <- x > 0 & y > 0
      out[valid] <- log2(y[valid] / x[valid])
      out
    }
    d$original_log2_bias <- ratio(d$input_freq, d$background_freq)
    d$sampled_log2_bias <- ratio(d$sampled_freq, d$background_freq)
    # Avoid Inf - Inf and 0/0 masquerading as a numerical improvement.
    d$bias_improvement <- abs(d$original_log2_bias) - abs(d$sampled_log2_bias)
    d$bias_improvement[!is.finite(d$bias_improvement)] <- NA_real_
    d$original_status <- ifelse(is.na(d$original_log2_bias), "Undefined",
                           ifelse(d$original_log2_bias > input$bias_cutoff, "Overrepresented",
                           ifelse(d$original_log2_bias < -input$bias_cutoff, "Underrepresented", "Within cutoff")))
    d$sampled_status <- ifelse(is.na(d$sampled_log2_bias), "Undefined",
                           ifelse(d$sampled_log2_bias > input$bias_cutoff, "Overrepresented",
                           ifelse(d$sampled_log2_bias < -input$bias_cutoff, "Underrepresented", "Within cutoff")))
    d
  })

  bias_display <- reactive({
    d <- bias_analysis()
    nmax <- input$top_n
    if (!is.null(nmax) && nmax > 0 && nrow(d) > nmax)
      d <- head(d, nmax)
    validate(need(nrow(d) > 0, "No groups to display."))
    # Stable cap shared by all three plots; infinite ratios remain in the table.
    finite <- c(d$original_log2_bias, d$sampled_log2_bias)
    finite <- abs(finite[is.finite(finite)])
    cap <- max(c(2, finite), na.rm=TRUE) + 0.5
    clamp <- function(x) {
      x[is.infinite(x) & x > 0] <- cap
      x[is.infinite(x) & x < 0] <- -cap
      x
    }
    d$original_display <- clamp(d$original_log2_bias)
    d$sampled_display <- clamp(d$sampled_log2_bias)
    d$display_cap <- cap
    d$display_label <- make.unique(paste0(d$name, " [", d$id, "]"))
    d
  })

  bias_colors <- c("Overrepresented"="#E68624", "Underrepresented"="#2479B5",
                   "Within cutoff"="#9B9B9B", "Undefined"="#D6D6D6")
  bias_tooltip <- function(d, idx) {
    paste0("Taxon: ", d$name[idx], "<br>ID: ", d$id[idx],
           "<br>Rank: ", d$rank[idx],
           "<br>Reference count: ", d$background[idx],
           "<br>Input count: ", d$input[idx],
           "<br>Sampled count: ", d$sampled[idx],
           "<br>Reference frequency: ", signif(d$background_freq[idx], 5),
           "<br>Input frequency: ", signif(d$input_freq[idx], 5),
           "<br>Sampled frequency: ", signif(d$sampled_freq[idx], 5),
           "<br>Original log2 bias: ", d$original_log2_bias[idx],
           "<br>Sampled log2 bias: ", d$sampled_log2_bias[idx],
           "<br>Bias improvement: ", d$bias_improvement[idx])
  }
  bias_long <- function(d) {
    idx <- rep(seq_len(nrow(d)), each=2)
    period <- rep(c("Original sequence DB", "Sampled sequences"), nrow(d))
    data.frame(taxon=d$display_label[idx],
               period=period,
               value=as.vector(t(as.matrix(d[,c("original_display","sampled_display")]))),
               status=as.vector(t(as.matrix(d[,c("original_status","sampled_status")]))),
               tooltip=paste0(bias_tooltip(d,idx), "<br>Dataset: ",period),
               stringsAsFactors=FALSE)
  }
  output$bias_ranked <- plotly::renderPlotly({
    d <- bias_display()
    long <- bias_long(d)
    long <- long[is.finite(long$value),,drop=FALSE]
    validate(need(nrow(long)>0, "No defined bias ratios to plot."))
    long$taxon <- factor(long$taxon, levels=rev(d$display_label))
    p <- ggplot(long, aes(x=taxon,y=value,fill=status,text=tooltip)) +
      geom_col(position=position_dodge(width=0.8),width=0.72) +
      facet_wrap(~period,nrow=1) + coord_flip() +
      geom_hline(yintercept=0,color="grey40") +
      scale_fill_manual(values=bias_colors,drop=FALSE) +
      labs(x=NULL,y="log2(dataset frequency / reference frequency)",fill="Bias classification") +
      theme_minimal(base_size=11) + theme(legend.position="bottom")
    plotly::ggplotly(p,tooltip="text")
  })
  output$bias_improvement <- plotly::renderPlotly({
    d <- bias_display()
    d <- d[is.finite(d$bias_improvement),,drop=FALSE]
    validate(need(nrow(d)>0,"No finite bias-improvement values. Check zero-frequency groups."))
    d <- d[order(d$bias_improvement),,drop=FALSE]
    d$display_label <- factor(d$display_label,levels=d$display_label)
    d$change <- ifelse(d$bias_improvement>0,"Improved",
                       ifelse(d$bias_improvement<0,"Deteriorated","Unchanged"))
    d$tooltip <- paste0(bias_tooltip(d,seq_len(nrow(d))),
                        "<br>Change: ",d$change)
    p <- ggplot(d,aes(x=display_label,y=bias_improvement,fill=change,text=tooltip)) +
      geom_col(width=0.75) + geom_hline(yintercept=0,color="grey40") +
      coord_flip() +
      scale_fill_manual(values=c("Improved"="#2479B5","Deteriorated"="#E68624","Unchanged"="#9B9B9B")) +
      labs(x=NULL,y="|Original bias| - |Sampled bias|",fill="Change") +
      theme_minimal(base_size=11) + theme(legend.position="bottom")
    plotly::ggplotly(p,tooltip="text")
  })
  output$bias_heatmap <- plotly::renderPlotly({
    d <- bias_display()
    long <- bias_long(d)
    long <- long[is.finite(long$value),,drop=FALSE]
    validate(need(nrow(long)>0,"No defined bias ratios to display."))
    long$taxon <- factor(long$taxon,levels=rev(d$display_label))
    long$period <- factor(long$period,levels=c("Original sequence DB","Sampled sequences"))
    # Diverging continuous scale: blue for negative, grey at zero, orange positive.
    p <- ggplot(long,aes(x=period,y=taxon,fill=value,text=tooltip)) +
      geom_tile(color="white",linewidth=0.35) +
      scale_fill_gradient2(low="#2479B5",mid="#9B9B9B",high="#E68624",
                           midpoint=0,limits=c(-d$display_cap[1],d$display_cap[1]),
                           name="log2 bias") +
      labs(x=NULL,y=NULL) + theme_minimal(base_size=11) +
      theme(panel.grid=element_blank(),axis.text.x=element_text(angle=15,hjust=1))
    plotly::ggplotly(p,tooltip="text")
  })
  output$download_bias <- downloadHandler(
    filename=function() paste0("TaxonSampling_bias_",Sys.Date(),".tsv"),
    content=function(file) utils::write.table(bias_analysis(),file,sep="\t",
                                               row.names=FALSE,quote=TRUE,na="NA")
  )

  # Retrieve parent/child rank relationships from the complete reference backbone.
  # A recursive CTE walks from each child-rank node to its nearest parent-rank
  # ancestor, including lineages absent from the sequence database.
  treemap_taxonomy <- reactive({
    req(state$db, state$counts)
    ranks <- c("phylum", "class", "order", "family", "genus", "species")
    root_rank <- input$tree_root_rank
    child_rank <- input$tree_child_rank
    validate(need(match(child_rank, ranks) > match(root_rank, ranks),
                  "Nested rank must be below the parent rank."))
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir=state$db, read_only=FALSE)
    on.exit(DBI::dbDisconnect(con, shutdown=TRUE), add=TRUE)
    # Do not stop recursion at a rank match: the independent diagnostic
    # established that doing so missed valid class-to-phylum relationships.
    # Traverse the complete ancestry, then select the closest matching parent.
    sql <- sprintf("WITH RECURSIVE ancestry AS (
      SELECT id AS child_id, id AS ancestor_id, parent AS next_parent,
             lower(trim(rank)) AS ancestor_rank, 0 AS depth
      FROM nodes WHERE lower(trim(rank)) = '%s'
      UNION ALL
      SELECT a.child_id, n.id, n.parent, lower(trim(n.rank)), a.depth + 1
      FROM ancestry a JOIN nodes n ON n.id = a.next_parent
      WHERE a.depth < 100 AND a.next_parent IS NOT NULL
        AND a.next_parent <> a.ancestor_id
    ),
    nearest_parent AS (
      SELECT child_id, ancestor_id AS root_id,
             ROW_NUMBER() OVER (PARTITION BY child_id ORDER BY depth) AS rn
      FROM ancestry WHERE depth > 0 AND ancestor_rank = '%s'
    ),
    nearest_same_rank AS (
      SELECT child_id, ancestor_id AS enclosing_child_id,
             ROW_NUMBER() OVER (PARTITION BY child_id ORDER BY depth) AS rn
      FROM ancestry WHERE depth > 0 AND ancestor_rank = '%s'
    )
    SELECT c.id AS child_id, c.name AS child_name,
           p.id AS root_id, p.name AS root_name,
           CASE WHEN e.enclosing_child_id IS NOT NULL
                     AND e.enclosing_child_id <> p.id
                THEN CAST(e.enclosing_child_id AS VARCHAR)
                ELSE NULL END AS enclosing_child_id,
           c.species_count AS child_species, c.tip_count AS child_tips,
           p.species_count AS root_species, p.tip_count AS root_tips
    FROM nodes c
    LEFT JOIN nearest_parent np ON np.child_id = c.id AND np.rn = 1
    LEFT JOIN nodes p ON p.id = np.root_id
    LEFT JOIN nearest_same_rank e ON e.child_id = c.id AND e.rn = 1
    WHERE lower(trim(c.rank)) = '%s'",
      child_rank, root_rank, child_rank, child_rank)
    mapped <- DBI::dbGetQuery(con, sql)
    if (anyDuplicated(as.character(mapped$child_id))) {
      stop('Treemap: duplicate child mapping; inspect ancestry query.')
    }
    mapped_count <- sum(!is.na(mapped$root_id))
    unassigned_count <- sum(is.na(mapped$root_id))
    message('Treemap: mapped ', mapped_count, ' ', child_rank,
            ' nodes to ', length(unique(mapped$root_id[!is.na(mapped$root_id)])),
            ' ', root_rank, ' nodes; unassigned: ', unassigned_count,
            '; total ', nrow(mapped), '.')
    # Explicit regression check for the synthetic phylum -> class case.
    # It is diagnostic, not a hard-coded assumption for other databases.
    if (root_rank == 'phylum' && child_rank == 'class' && nrow(mapped) == 249L &&
        (mapped_count != 184L || unassigned_count != 65L)) {
      warning(sprintf('Synthetic treemap regression: expected 184 mapped / 65 unassigned classes; got %d / %d.',
                      mapped_count, unassigned_count))
    }
    roots <- DBI::dbGetQuery(con, sprintf(
      "SELECT id AS root_id, name AS root_name, species_count AS root_species, tip_count AS root_tips FROM nodes WHERE lower(trim(rank)) = '%s'", root_rank))
    list(roots=roots, children=mapped)
  })

  # get_counts() intentionally prunes the scaffold to observed nodes.
  # A class with no observations is absent from that object, even when its
  # parent phylum has observations. Therefore counts for treemaps MUST be
  # propagated on the full DuckDB hierarchy, not looked up in comparison_data().
  treemap_observed_counts <- reactive({
    req(state$db, state$counts)
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir=state$db, read_only=FALSE)
    on.exit(DBI::dbDisconnect(con, shutdown=TRUE), add=TRUE)
    ranks <- c('phylum','class','order','family','genus','species')
    parent_rank <- input$tree_root_rank
    child_rank <- input$tree_child_rank
    validate(need(parent_rank %in% ranks && child_rank %in% ranks &&
                  match(child_rank,ranks) > match(parent_rank,ranks),
                  'Select a nested rank below the parent rank.'))

    # One row per sequence; repeated taxon IDs must remain repeated.
    ids <- as.character(state$counts$ids_df$node_id)
    ids <- ids[!is.na(ids) & nzchar(ids)]
    sampled <- if (!is.null(state$result)) as.character(state$result$outputNodes) else character()
    sampled <- sampled[!is.na(sampled) & nzchar(sampled)]
    entries <- rbind(
      data.frame(node_id=ids, dataset='input', stringsAsFactors=FALSE),
      data.frame(node_id=sampled, dataset='sampled', stringsAsFactors=FALSE)
    )
    if (!nrow(entries)) return(data.frame(node_id=character(),dataset=character(),n=numeric()))
    entries <- as.data.frame(table(entries$node_id, entries$dataset), stringsAsFactors=FALSE)
    names(entries) <- c('node_id','dataset','n')
    entries <- entries[entries$n > 0, , drop=FALSE]
    entries$node_id <- as.character(entries$node_id)
    entries$dataset <- as.character(entries$dataset)
    entries$n <- as.numeric(entries$n)
    DBI::dbWriteTable(con, 'ts_treemap_input', entries, temporary=TRUE, overwrite=TRUE)
    # Traverse each observed node upwards. Include its own rank and every
    # ancestor, but stop at the root; aggregate by node ID and dataset.
    sql <- sprintf("WITH RECURSIVE ancestors AS (
      SELECT i.node_id AS id, i.dataset, i.n, n.parent, n.rank, 0 AS depth
      FROM ts_treemap_input i JOIN nodes n ON n.id=i.node_id
      UNION ALL
      SELECT p.id, a.dataset, a.n, p.parent, p.rank, a.depth+1
      FROM ancestors a JOIN nodes p ON p.id=a.parent
      WHERE a.depth < 100 AND a.parent IS NOT NULL AND a.parent <> a.id
    )
    SELECT id AS node_id, dataset, SUM(n) AS n
    FROM ancestors WHERE lower(rank) IN ('%s','%s')
    GROUP BY id, dataset", parent_rank, child_rank)
    DBI::dbGetQuery(con, sql)
  })

  treemap_data <- reactive({
    req(state$counts)
    taxonomy <- treemap_taxonomy()
    roots <- taxonomy$roots
    children <- taxonomy$children
    validate(need(nrow(roots) > 0, 'No parent-rank groups in the reference database.'))
    observed <- treemap_observed_counts()
    unit <- state$counts$metadata$background_unit
    bg_col <- if (identical(unit, 'tips')) 'tips' else 'species'
    lookup <- function(ids, dataset) {
      rows <- observed[observed$dataset == dataset, , drop=FALSE]
      v <- rows$n[match(as.character(ids), as.character(rows$node_id))]
      v[is.na(v)] <- 0
      as.numeric(v)
    }
    palette <- c('#636EFA','#EF553B','#00CC96','#AB63FA',
                 '#FFA15A','#19D3F3','#FF6692','#B6E880','#FF97FF','#FECB52')
    root_keys <- sort(unique(as.character(roots$root_id)))
    root_colors <- setNames(rep(palette,length.out=length(root_keys)),root_keys)
    result <- list()
    for (kind in c('background','input','sampled')) {
      rv <- if (kind == 'background') as.numeric(roots[[paste0('root_',bg_col)]]) else
        lookup(roots$root_id, kind)
      cv <- if (kind == 'background') as.numeric(children[[paste0('child_',bg_col)]]) else
        lookup(children$child_id, kind)
      # Include all positive classes, including the ones with no phylum ancestor.
      # Same-rank nesting is represented by an actual parent-child edge.
      active_root <- which(is.finite(rv) & rv > 0)
      active_child <- which(is.finite(cv) & cv > 0)
      root_ids <- as.character(roots$root_id[active_root])
      child_ids <- as.character(children$child_id[active_child])
      assigned_root <- as.character(children$root_id[active_child])
      enclosing <- as.character(children$enclosing_child_id[active_child])
      enclosing[is.na(enclosing)] <- ''
      # A same-rank ancestor is only a valid direct parent if it belongs to
      # the same root group (including both being unassigned).
      ancestor_index <- match(enclosing, as.character(children$child_id))
      ancestor_root <- as.character(children$root_id[ancestor_index])
      valid_enclosing <- nzchar(enclosing) & enclosing %in% child_ids &
        ((is.na(assigned_root) & is.na(ancestor_root)) |
         (!is.na(assigned_root) & !is.na(ancestor_root) & assigned_root == ancestor_root))
      child_parent <- ifelse(valid_enclosing, paste0('child:',enclosing),
                             ifelse(is.na(assigned_root), 'root:unassigned',
                                    paste0('root:',assigned_root)))
      # Ensure the parent exists; otherwise report rather than silently drop.
      missing_roots <- setdiff(unique(assigned_root[!is.na(assigned_root)]),root_ids)
      if (length(missing_roots)) {
        stop(sprintf('Treemap (%s): %d positive child groups have absent parent groups: %s',
                     kind,length(missing_roots),paste(head(missing_roots,8),collapse=', ')))
      }
      has_unassigned <- any(is.na(assigned_root))
      ids <- c(paste0('root:',root_ids),if (has_unassigned) 'root:unassigned',
               paste0('child:',child_ids))
      labels <- c(as.character(roots$root_name[active_root]),
                  if (has_unassigned) paste('Unassigned to',input$tree_root_rank),
                  as.character(children$child_name[active_child]))
      # For unassigned taxa, derive the synthetic root from the sum of
      # top-level (non-nested) child totals; never invent a phylum assignment.
      unassigned_total <- sum(cv[active_child][is.na(assigned_root) &
        child_parent == 'root:unassigned'],na.rm=TRUE)
      values <- c(rv[active_root],if (has_unassigned) unassigned_total,
                  cv[active_child])
      parents <- c(rep('',length(active_root)),if (has_unassigned) '',child_parent)
      colors <- c(unname(root_colors[root_ids]),if (has_unassigned) '#888888',
                  ifelse(is.na(assigned_root),'#888888',
                         unname(root_colors[assigned_root])))
      # Plotly 'total' requires the immediate children to fit inside each parent.
      # Do not omit children or rescale data to hide an inconsistency.
      child_sums <- tapply(values[nzchar(parents)],parents[nzchar(parents)],sum)
      node_values <- setNames(values,ids)
      excess <- child_sums - node_values[names(child_sums)]
      bad <- names(excess)[!is.finite(excess) |
        excess > 1e-8*pmax(1,node_values[names(excess)])]
      if (length(bad)) {
        detail <- paste(head(sprintf('%s: children=%g, parent=%g',
                         bad,child_sums[bad],node_values[bad]),8),collapse='; ')
        stop(sprintf('Treemap (%s): overlapping or inconsistent counts: %s',kind,detail))
      }
      # Only selected parent-rank nodes receive an explicit residual.
      # Child-rank nodes already represent the complete descendant count;
      # a residual inside a terminal child would misleadingly label its
      # unshown orders/families/species as "unclassified".
      # Nested same-rank child nodes remain nested: Plotly draws the enclosing
      # child's undisplayed remainder as part of that child's own rectangle.
      parent_ids <- paste0('root:', root_ids)
      parent_remaining <- setNames(rv[active_root], parent_ids)
      root_child_sums <- child_sums[intersect(names(child_sums), parent_ids)]
      if (length(root_child_sums)) {
        parent_remaining[names(root_child_sums)] <- pmax(
          0, parent_remaining[names(root_child_sums)] - root_child_sums
        )
      }
      residual <- names(parent_remaining)[
        is.finite(parent_remaining) & parent_remaining > 1e-8
      ]
      ids <- c(ids, paste0('other:', residual))
      labels <- c(labels, rep('Other / unclassified', length(residual)))
      parents <- c(parents, residual)
      values <- c(values, unname(parent_remaining[residual]))
      colors <- c(colors, rep('#A3A3A3', length(residual)))
      result[[kind]] <- data.frame(ids=ids,labels=labels,parents=parents,
                                   values=values,colors=colors,
                                   stringsAsFactors=FALSE)
    }
    result
  })

  draw_treemap <- function(kind, title) {
    d <- treemap_data()[[kind]]
    validate(need(!is.null(d) && nrow(d) > 0,
                  paste("No", title, "counts available at the selected ranks.")))
    plotly::plot_ly(type="treemap", ids=d$ids, labels=d$labels,
                    parents=d$parents, values=d$values,
                    marker=list(colors=d$colors),
                    branchvalues="total", textinfo="label+percent parent",
                    hovertemplate=paste0("%{label}<br>Count: %{value:,}",
                                         "<br>Share of parent: %{percentParent:.2%}<extra></extra>")) |>
      plotly::layout(title=list(text=title), margin=list(l=5,r=5,b=5,t=45))
  }
  output$tree_reference <- plotly::renderPlotly({ draw_treemap("background", "Full reference database") })
  output$tree_input <- plotly::renderPlotly({ draw_treemap("input", "Full sequence database") })
  output$tree_sampled <- plotly::renderPlotly({ req(state$result); draw_treemap("sampled", "Sampled sequences") })

  output$download_counts <- downloadHandler(filename=function() "taxonsampling_counts.rds", content=function(file) {req(state$counts); saveRDS(state$counts, file)})
  output$download_result <- downloadHandler(filename=function() "taxonsampling_result.rds", content=function(file) {req(state$result); saveRDS(state$result, file)})
  output$download_selected <- downloadHandler(filename=function() "taxonsampling_selected.tsv", content=function(file) {req(state$result); write.table(state$result$outputMapping, file, sep="\t", quote=FALSE, row.names=FALSE)})
  output$download_summary <- downloadHandler(filename=function() "taxonsampling_comparison.tsv", content=function(file) {write.table(filtered_comparison(), file, sep="\t", quote=FALSE, row.names=FALSE)})
}

shinyApp(ui, server)
