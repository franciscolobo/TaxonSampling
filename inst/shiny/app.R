library(shiny)
library(TaxonSampling)
library(DT)
library(ggplot2)
library(plotly)
library(shinyFiles)

# Paths to local reference files are entered on the host running Shiny.
ui <- fluidPage(
  titlePanel("TaxonSampling | phylogeny-aware sampling"),
  sidebarLayout(
    sidebarPanel(width = 4,
      tags$img(src = "TS_logo.png", style = "max-width:180px;max-height:100px;"),
      tabsetPanel(id = "workflow",
        tabPanel("1. Backbone",
          selectInput("taxonomy", "Reference", c("NCBI Taxonomy"="ncbi", "Open Tree Taxonomy"="ott", "Open Tree synthetic phylogeny"="synthetic", "Custom Newick tree"="custom")),
          textInput("taxonomy_path", "Extracted taxonomy directory / Newick file (server path)"),
          shinyFilesButton("browse_taxonomy_file", "Browse file", "Select Newick or taxonomy file", multiple=FALSE),
          shinyDirButton("browse_taxonomy_dir", "Browse directory", "Select extracted taxonomy directory"),
          conditionalPanel("input.taxonomy == 'synthetic'",
            textInput("ott_path", "OTT directory for synthetic species annotations (optional)"),
            shinyDirButton("browse_ott_dir", "Browse directory", "Select OTT directory")),
          textInput("output_dir", "Database output directory", "taxonomy_db"),
          shinyDirButton("browse_output_dir", "Browse directory", "Select database output directory"),
          textInput("db_name", "Database filename (optional; defaults to backbone name)", ""),
          checkboxInput("include_suppressed", "Include suppressed OTT taxa", TRUE),
          checkboxInput("overwrite", "Overwrite existing database", FALSE),
          actionButton("build", "Build reference database"),
          tags$hr(),
          textInput("db_path", "Existing DuckDB file (server path)", "taxonomy_db/ncbi.duckdb"),
          shinyFilesButton("browse_db_file", "Browse database", "Select DuckDB database", multiple=FALSE),
          actionButton("inspect", "Inspect database")
        ),
        tabPanel("2. Input taxa",
          fileInput("taxa_file", "Taxon / sequence ID TSV (headerless, two columns)", accept=c(".tsv", ".txt", ".csv")),
          helpText("One row per taxon–sequence pair. IDs are strings; retain leading zeros. Use the same taxon with multiple sequence IDs if needed."),
          selectInput("id_type", "Identifier namespace", c("NCBI"="ncbi", "OTT"="ott", "Native backbone ID"="native")),
          selectInput("count_mode", "Observation unit", c("Distinct input taxon IDs"="input_ids", "Distinct mapped backbone nodes"="backbone_nodes")),
          selectInput("background_unit", "Background unit", c("Species"="species", "Terminal tips"="tips")),
          selectInput("missing_ids", "Unmapped identifiers", c("Warn"="warn", "Drop"="drop", "Error"="error")),
          selectInput("duplicate_ids", "Repeated identical taxon/sequence pairs", c("Keep first"="first", "Error"="error")),
          actionButton("count", "Prepare dataset counts"),
          downloadButton("download_counts", "Download count object (.rds)")
        ),
        tabPanel("3. Sampling",
          textInput("root", "Root node ID (native backbone identifier)", ""),
          helpText("Use an ID from the observed nodes table. OTT IDs in synthetic trees may use the ott: prefix."),
          numericInput("m", "Number of observations", 100, min=1, step=1),
          selectInput("method", "Method", c("Diversity"="diversity", "Balanced (experimental)"="balanced")),
          selectInput("randomize", "Randomization", c("After first round"="after_first_round", "Yes"="yes", "No"="no")),
          selectInput("sampling", "Sampling mode", c("Known species"="known_species", "Agnostic"="agnostic")),
          checkboxInput("replacement", "Sample with replacement", FALSE),
          textAreaInput("ignore", "Exclude clade IDs (comma / newline separated)", "", rows=2),
          textAreaInput("require", "Require clade IDs (comma / newline separated)", "", rows=2),
          numericInput("seed", "Random seed", 1, min=0, step=1),
          actionButton("sample", "Run TaxonSampling"),
          downloadButton("download_selected", "Download sampled taxa (.tsv)"),
          downloadButton("download_result", "Download full result (.rds)")
        ),
        tabPanel("4. Visualization",
          selectInput("level", "Taxonomic rank", c("Phylum"="phylum", "Class"="class", "Order"="order", "Family"="family", "Genus"="genus", "Species"="species", "All observed nodes"="all")),
          numericInput("top_n", "Maximum plotted groups (0 = all)", 40, min=0, max=100000),
          numericInput("bias_cutoff", "Bias cutoff: absolute log2 fold change", 1, min=0, step=0.25),
          helpText("Orange = overrepresented (Y versus X); blue = underrepresented; grey = within cutoff. Cutoff 1 = twofold. This is an effect-size threshold, not a statistical significance test."),
          checkboxInput("labels", "Show taxon labels", FALSE),
          downloadButton("download_summary", "Download comparison table (.tsv)")
        )
      )
    ),
    mainPanel(width=8,
      verbatimTextOutput("status"),
      tabsetPanel(
        tabPanel("Reference", DTOutput("reference_table")),
        tabPanel("Observed nodes", DTOutput("nodes_table")),
        tabPanel("Mappings", DTOutput("mapping_table")),
        tabPanel("Selected taxa", DTOutput("selected_table")),
        tabPanel("Representation",
          h4("1. Full reference vs full input database"), plotlyOutput("comparison_reference_input", height="490px"),
          h4("2. Full reference vs sampled dataset"), plotlyOutput("comparison_reference_sample", height="490px"),
          h4("3. Full input database vs sampled dataset"), plotlyOutput("comparison_input_sample", height="490px"),
          DTOutput("summary_table")),
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
  state <- reactiveValues(db=NULL, counts=NULL, result=NULL, message="Choose an existing reference database or build one.")
  report_error <- function(e) { state$message <- paste("Error:", conditionMessage(e)); showNotification(state$message, type="error", duration=NULL) }
  parse_ids <- function(s) { x <- trimws(unlist(strsplit(s, "[,;\\n\\r\\t ]+"))); unique(x[nzchar(x)]) }
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
      updateTextInput(session, "db_path", value=state$db)
      state$counts <- NULL; state$result <- NULL
      state$message <- paste("Reference database created:", state$db)
    }, error=report_error)
  })
  observeEvent(input$inspect, {
    tryCatch({
      path <- normalizePath(trimws(input$db_path), mustWork=TRUE)
      con <- DBI::dbConnect(duckdb::duckdb(), dbdir=path, read_only=TRUE)
      on.exit(DBI::dbDisconnect(con, shutdown=TRUE))
      if (!all(c("nodes", "metadata") %in% DBI::dbListTables(con))) stop("Not a TaxonSampling reference database")
      state$db <- path; state$counts <- NULL; state$result <- NULL
      state$message <- paste("Reference database loaded:", path)
    }, error=report_error)
  })
  reference <- reactive({
    req(state$db)
    con <- DBI::dbConnect(duckdb::duckdb(), dbdir=state$db, read_only=TRUE)
    on.exit(DBI::dbDisconnect(con, shutdown=TRUE))
    DBI::dbGetQuery(con, "SELECT rank, COUNT(*) AS nodes, SUM(CASE WHEN node_type = 'tip' THEN 1 ELSE 0 END) AS tips FROM nodes GROUP BY rank ORDER BY nodes DESC LIMIT 100")
  })
  output$reference_table <- renderDT(datatable(reference(), options=list(pageLength=15, scrollX=TRUE)))
  observeEvent(input$count, {
    tryCatch({
      req(state$db, input$taxa_file)
      # File upload path is transient; get_counts reads it immediately.
      counts <- TaxonSampling::get_counts(
        taxonomy_db=state$db, ids_file=input$taxa_file$datapath,
        id_type=input$id_type, count_mode=input$count_mode,
        background_unit=input$background_unit, missing_ids=input$missing_ids,
        duplicate_ids=input$duplicate_ids
      )
      state$counts <- counts; state$result <- NULL
      roots <- counts$nodes$id[is.na(counts$nodes$parent) | !(counts$nodes$parent %in% counts$nodes$id) | counts$nodes$parent == counts$nodes$id]
      if (length(roots)) updateTextInput(session, "root", value=as.character(roots[1]))
      state$message <- sprintf("Counts prepared: %s mapped taxon/sequence pairs; %s observed nodes; %s unmapped IDs.", nrow(counts$ids_df), nrow(counts$nodes), length(counts$unmapped_ids))
    }, error=report_error)
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
  output$download_counts <- downloadHandler(filename=function() "taxonsampling_counts.rds", content=function(file) {req(state$counts); saveRDS(state$counts, file)})
  output$download_result <- downloadHandler(filename=function() "taxonsampling_result.rds", content=function(file) {req(state$result); saveRDS(state$result, file)})
  output$download_selected <- downloadHandler(filename=function() "taxonsampling_selected.tsv", content=function(file) {req(state$result); write.table(state$result$outputMapping, file, sep="\t", quote=FALSE, row.names=FALSE)})
  output$download_summary <- downloadHandler(filename=function() "taxonsampling_comparison.tsv", content=function(file) {write.table(filtered_comparison(), file, sep="\t", quote=FALSE, row.names=FALSE)})
}

shinyApp(ui, server)
