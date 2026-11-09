library(shiny)
library(shinyFiles)
library(fs)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(TaxonSampling)
#devtools::load_all("../TaxonSampling/")
library(tidyr)
library(magrittr)
library(dplyr)
library(plotly)
library(DT)


shinyServer(function(input, output, session) {

  #Starts volumes and files
  volumes <- c(NCBI_Tax = "/Users/chico/projects/TaxonSampling/data/NCBI_Taxonomy_data/", TS = "/Users/chico/projects/TaxonSampling/data/TS_data/metadata/", Output = "/Users/chico/antigo_F/lala/", getVolumes()())

  #User-provided tabular file linking NCBI Taxonomy IDs to sequence IDs
  shinyFileChoose(input, "user_database_file_path", roots = volumes, session = session)
  #Path to NCBI Taxonomy dump dir
  shinyDirChoose(input, "NCBI_Taxonomy_directory", roots = volumes, session = session, restrictions = system.file(package = "base"), allowDirCreate = FALSE)
  #Path to write count database
  shinyDirChoose(input, "TS_DB_database_path", roots = volumes, session = session, restrictions = system.file(package = "base"), allowDirCreate = TRUE)
  #Path to previously saved count database
  shinyDirChoose(input, "TS_DB_database_path_previous_run", roots = volumes, session = session, restrictions = system.file(package = "base"), allowDirCreate = TRUE)
  #Path to input
  shinyDirChoose(input, "TS_DB_database_input_path", roots = volumes, session = session, restrictions = system.file(package = "base"), allowDirCreate = TRUE)
  
  shinyDirChoose(input, "TS_outptut_path", roots = volumes, session = session, restrictions = system.file(package = "base"), allowDirCreate = TRUE)

  shinyDirChoose(input, "TS_DB_input_path", roots = volumes, session = session, restrictions = system.file(package = "base"), allowDirCreate = TRUE)
  
  shinyDirChoose(input, "TS_run_input_path", roots = volumes, session = session, restrictions = system.file(package = "base"), allowDirCreate = TRUE)
  
  
  
  ## print to console to see how the value of the shinyFiles 
  ## button changes after clicking and selection

  #uptdate UI output tab
  observeEvent(input$task, {
    updateTabsetPanel(inputId = "outtabs", selected = input$task)
  })
  
  #uptdate UI input parameters
  observeEvent(input$task, {
    updateTabsetPanel(inputId = "params", selected = input$task)
  }) 
  
  
  #start making a countDB database, or using a previously computed one
  
  #create countDB
  count_db <- eventReactive(input$start_analysis_DB_build, {
    isolate({
      ids_file_path <- parseFilePaths(volumes, input$user_database_file_path)
#      ids_file_path
#      class(ids_file_path)
      output$ids_file_path <- renderPrint(ids_file_path)
      NCBI_Tax_input_dir_path <- parseDirPath(volumes, input$NCBI_Taxonomy_directory)
      #      get_counts(taxonomy_path = input_dir_path, ids_file = "/Users/chico/projects/TaxonSampling/data/TS_data/metadata/TaxID2SeqID.txt", start_from_species = TRUE)
      count_db <- get_counts(taxonomy_path = NCBI_Tax_input_dir_path, ids_file = ids_file_path$datapath, start_from_species = TRUE)
      out_file_path <- paste0(parseDirPath(volumes, input$TS_DB_database_path), "/", "count_db.RData")
      save(count_db, file = out_file_path)
    })
  })
  
  output$create_countsDB <- renderPrint({
    db <- count_db()
    print("Done!")
  })

  output$create_TSDB <- renderPrint({
    TS <- run_TaxonSampling()
    print("Done!")
  })
  

  run_TaxonSampling <- eventReactive(input$start_analysis_TS, {
    isolate({
      input_file_path <- paste0(parseDirPath(volumes, input$TS_DB_database_input_path), "/", "count_db.RData")
#      input_file_path
      load(input_file_path)
      TS_output <- run_TS(taxlist = count_db, taxon = input$root_taxon, m = input$m, method = input$method, randomize = input$randomize, sampling = input$sampling)
      output_file_path <- paste0(parseDirPath(volumes, input$TS_outptut_path), "/", "TS_output.RData")
      save(TS_output, file = output_file_path)
    })
  })
  
  
  ## print directory and file paths to browser
  
  output$TaxID2SeqID_file_path <- renderPrint({
    if (is.integer(input$user_database_file_path)) {
      cat("No files have been selected")
    } else {
      file_path <- parseFilePaths(volumes, input$user_database_file_path)
      as.character(file_path$datapath)
    }
  })
  
  output$NCBI_Tax_path <- renderPrint({
    if (is.integer(input$NCBI_Taxonomy_directory)) {
      cat("No directory selected")
    } else {
      parseDirPath(volumes, input$NCBI_Taxonomy_directory)
    }
  })
  
  output$saveDir_path <- renderPrint({
    if(is.integer(input$TS_DB_database_path)) {
      cat("No directory selected")
    } else {
      parseDirPath(volumes, input$TS_DB_database_path)
    }
  })

  output$saveDir_TS_path <- renderPrint({
    if(is.integer(input$TS_outptut_path)) {
      cat("No directory selected")
    } else {
      parseDirPath(volumes, input$TS_outptut_path)
    }
  })

  output$inputDir_TS_path <- renderPrint({
    if(is.integer(input$TS_DB_database_input_path)) {
      cat("No directory selected")
    } else {
      parseDirPath(volumes, input$TS_DB_database_input_path)
    }
  })
  
  output$saveDir_path <- renderPrint({
    if(is.integer(input$TS_DB_database_path)) {
      cat("No directory selected")
    } else {
      parseDirPath(volumes, input$TS_DB_database_path)
    }
  })

  output$TS_DB_input_path <- renderPrint({
    if(is.integer(input$TS_DB_input_path)) {
      cat("No directory selected")
    } else {
      parseDirPath(volumes, input$TS_DB_input_path)
    }
  })
  
  output$TS_run_input_path <- renderPrint({
    if(is.integer(input$TS_run_input_path)) {
      cat("No directory selected")
    } else {
      parseDirPath(volumes, input$TS_run_input_path)
    }
  })

  output$DB_summary <- renderDataTable( 
    df <- get_DB_summary(),
    options = list(scrollX = TRUE),
    rownames = FALSE
  )

  get_DB_summary <- eventReactive(input$start_analysis_DB_plot, {
    input_db_file_path <- paste0(parseDirPath(volumes, input$TS_DB_database_path_previous_run), "/", "count_db.RData")
    #load "count_db" object in the environment
    load(input_db_file_path)
  
    #define root taxon to filter data
    root_taxon <- input$root_taxon
  
    #get summary (maybe move this function to TS?)
    summary <- get_summary_count_db(count_db, root_taxon)
    
    #start building plot df
    #summarize by user-defined level
    df <- summary[summary$level == input$taxonomic_level,]
    df$species_freq <- df$species_count/sum(df$species_count)
    #frequency of spp in user-defined db
    df$seq_freq <- df$seq_count/sum(df$seq_count)
    #log2fold change of frequencies of databasse spp and known spp
    df$log2_fold_change <- round(log2(df$seq_freq) - log2(df$species_freq), digits = 2)
    df$species_freq <- formatC((df$species_count/sum(df$species_count)), format = "E", digits = 2)
    #frequency of spp in user-defined db
    df$seq_freq <- formatC((df$seq_count/sum(df$seq_count)), format = "E", digits = 2)
    df
#    datatable(df) %>%
#      formatSignif(columns = c("species_freq", "seq_freq"),
#                   digits = 2)
    
  })
  
  get_scatterplot_DB <- eventReactive(input$start_analysis_DB_plot, {
    
    #path to previously build database
    input_db_file_path <- paste0(parseDirPath(volumes, input$TS_DB_database_path_previous_run), "/", "count_db.RData")
    
    #load "count_db" object in the environment
    load(input_db_file_path)
    
    root_taxon <- input$root_taxon
    
    #get summary (maybe move this function to TS?)
    summary <- get_summary_count_db(count_db, root_taxon)

    #start building plot df
    #summarize by user-defined level
    plot_df <- summary[summary$level == input$taxonomic_level,]
    #frequency of knonw ssp
    plot_df$species_freq <- plot_df$species_count/sum(plot_df$species_count)
    #frequency of spp in user-defined db
    plot_df$seq_freq <- plot_df$seq_count/sum(plot_df$seq_count)
    #log2fold change of frequencies of databasse spp and known spp
    plot_df$log2_fold_change <- round(log2(plot_df$seq_freq) - log2(plot_df$species_freq), digits = 2)
    #defining colors
    plot_df$color <- "gray"
    #user-defined upper bound for log2fc
    plot_df$color[plot_df$log2_fold_change > input$upper_bound] <- "orange"
    #user-defined lower bound for log2fc
    plot_df$color[plot_df$log2_fold_change < input$lower_bound] <- "blue"
    #plot db to compute regression (values within upper and lower bound cutoffs)
    plot_df_regression <- plot_df[(plot_df$log2_fold_change > input$lower_bound & plot_df$log2_fold_change < input$upper_bound), ]

    #plotting
    plot <- plot_df %>%
      
      ggplot(aes(x=species_freq, y=seq_freq)) +
      
      #adding correlation for linear regression
      stat_cor(data=plot_df_regression, aes(label = as.character(after_stat(r.label))),
        method = "pearson",
        label.x.npc = "left",
        label.y.npc = "top"
      ) +
      
      geom_smooth(method=lm , color="gray", se=TRUE) + 
      
      geom_point(aes(color=color,
                  text=paste("Name: ", name, "<br>", 
                             "Known species count: ", species_count, "<br>",
                             "Known species frequency: ", species_freq, "<br>",
                             "DB species count: ", seq_count, "<br>",
                             "DB species frequency: ", seq_freq, "<br>",
                             "log2 fold change: ", log2_fold_change,
                             sep = "")))+
  
      
      scale_color_identity() +
      
#      scale_color_gradient2(low = "#0096FF", mid = "#F7F7F7", high ="#F1A340", 
#                           midpoint = 0, space = "Lab", guide = "colourbar") +
      
      scale_size_continuous(name = "Number of\nknown species",
                            breaks = c(1, 100, 1000, 10000, 25000, 50000),
                            limits = c(0, max(plot_df$species_count))) + 
      
      scale_x_log10(labels = scales::scientific) +
      scale_y_log10(labels = scales::scientific) +
      
      xlab("Known species") +
      ylab("Species in database") +
      
      theme_light() +
      theme(
        axis.text=element_text(size=rel(1)),
        axis.title=element_text(size=rel(1))
      )

    plot <- plot + geom_smooth(data=plot_df_regression, aes(x=species_freq, y=seq_freq), method=lm , color="red", se=TRUE)
    
    ggplotly(plot, tooltip = c("text"))
  })

  
  output$TS_summary <- renderDataTable( 
    df <- get_TS_summary(),
    options = list(scrollX = TRUE),
    rownames = FALSE
  )
  
  get_TS_summary <- eventReactive(input$start_analysis_TS_plot, {
    input_TS_file_path <- paste0(parseDirPath(volumes, input$TS_run_input_path), "/", "TS_output.RData")
    load(input_TS_file_path)
    root_taxon <- TS_output$ts.params$taxon
    
    summary_ts <- get_summary_df(TS_output)
    
    #start building plot df
    #summarize by user-defined level
    plot_df <- summary_ts[summary_ts$level == input$taxonomic_level_TS,]
    plot_df$sample_freq <- formatC(plot_df$sample_count/sum(plot_df$sample_count), format = "E", digits = 2)
    #frequency of knonw ssp
    plot_df$species_freq <- formatC(plot_df$species_count/sum(plot_df$species_count), format = "E", digits = 2)
    #frequency of spp in user-defined db
    plot_df$seq_freq <- formatC(plot_df$seq_count/sum(plot_df$seq_count), format = "E", digits = 2)
    #log2fold change of frequencies of databasse spp and known spp
    plot_df$log2FC_DB_vs_Spp <- round(log2(plot_df$seq_freq) - log2(plot_df$species_freq), digits = 2)
    plot_df$log2FC_TS_vs_Spp <- round(log2(plot_df$sample_freq) - log2(plot_df$species_freq), digits = 2)
    #defining colors
#    plot_df$color <- "gray"
    #user-defined upper bound for log2fc
#    plot_df$color[plot_df$log2_fold_change > input$upper_bound_TS] <- "orange"
    #user-defined lower bound for log2fc
#    plot_df$color[plot_df$log2_fold_change < input$lower_bound_TS] <- "blue"
    #plot db to compute regression (values within upper and lower bound cutoffs)
#    plot_df_regression <- plot_df[(plot_df$log2_fold_change > input$lower_bound_TS & plot_df$log2_fold_change < input$upper_bound_TS), ]
    plot_df
    })  
  
  

get_scatterplot_DB_TS <- eventReactive(input$start_analysis_TS_plot, {
  
  #path to TS output (needed to get root taxon to start displaying data)
  input_TS_file_path <- paste0(parseDirPath(volumes, input$TS_run_input_path), "/", "TS_output.RData")
  
  load(input_TS_file_path)
  
  root_taxon <- TS_output$ts.params$taxon
  
  print(root_taxon)
  
  rm(TS_output)
  
  #path to previously build database
  input_db_file_path <- paste0(parseDirPath(volumes, input$TS_DB_input_path), "/", "count_db.RData")
  
  #load "count_db" object in the environment
  load(input_db_file_path)
  
  #get summary (maybe move this function to TS?)
  summary <- get_summary_count_db(count_db, root_taxon)
  
  #start building plot df
  #summarize by user-defined level
  plot_df <- summary[summary$level == input$taxonomic_level_TS,]
  #frequency of knonw ssp
  plot_df$species_freq <- plot_df$species_count/sum(plot_df$species_count)
  #frequency of spp in user-defined db
  plot_df$seq_freq <- plot_df$seq_count/sum(plot_df$seq_count)
  #log2fold change of frequencies of databasse spp and known spp
  plot_df$log2_fold_change <- round(log2(plot_df$seq_freq) - log2(plot_df$species_freq), digits = 2)
  #defining colors
  plot_df$color <- "gray"
  #user-defined upper bound for log2fc
  plot_df$color[plot_df$log2_fold_change > input$upper_bound_TS] <- "orange"
  #user-defined lower bound for log2fc
  plot_df$color[plot_df$log2_fold_change < input$lower_bound_TS] <- "blue"
  #plot db to compute regression (values within upper and lower bound cutoffs)
  plot_df_regression <- plot_df[(plot_df$log2_fold_change > input$lower_bound_TS & plot_df$log2_fold_change < input$upper_bound_TS), ]
  
  #plotting
  plot <- plot_df %>%
    
    ggplot(aes(x=species_freq, y=seq_freq)) +
    
    #adding correlation for linear regression
    stat_cor(data=plot_df_regression, aes(label = as.character(after_stat(r.label))),
             method = "pearsonPearson",
             label.x.npc = "left",
             label.y.npc = "top"
    ) +
    
    #      geom_smooth(method=lm , color="gray", se=TRUE) + 
    
    geom_point(aes(name=name, log2_fold_change=log2_fold_change, size=species_count, seq_count=seq_count, color=color)) + 
    #                 text= paste("Name: ", Labels, "<br>", 
    #                             "Known species count: ", species_count, "<br>",
    #                             "Known species frequency: ", species_freq, "<br>",
    #                             "DB species frequency: ", seq_freq,
    #                             "log2 fold change: ", `log2 fold change`,
    #                             sep = ""))+
    
    
    scale_color_identity() +
    
    #      scale_color_gradient2(low = "#0096FF", mid = "#F7F7F7", high ="#F1A340", 
    #                           midpoint = 0, space = "Lab", guide = "colourbar") +
    
    scale_size_continuous(name = "Number of\nknown species",
                          breaks = c(1, 100, 1000, 10000, 25000, 50000),
                          limits = c(0, max(plot_df$species_count))) + 
    
    scale_x_log10(labels = scales::scientific) +
    scale_y_log10(labels = scales::scientific) +
    
    xlab("Known species") +
    ylab("Species in database") +
    
    theme_light() +
    theme(
      axis.text=element_text(size=rel(1)),
      axis.title=element_text(size=rel(1))
    )
  
  plot <- plot + geom_smooth(data=plot_df_regression, aes(x=species_freq, y=seq_freq), method=lm , color="red", se=TRUE)
  
  ggplotly(plot, tooltip = c("name", "log2_fold_change", "species_count", "seq_count"))
})



get_scatterplot_TS <- eventReactive(input$start_analysis_TS_plot, {

  input_TS_file_path <- paste0(parseDirPath(volumes, input$TS_run_input_path), "/", "TS_output.RData")
  
  load(input_TS_file_path)
  
  summary <- get_summary_df(TS_output)

  plot_df <- summary[summary$level == input$taxonomic_level_TS,]
  
  #})
  
  #get_scatterplot2 <- eventReactive(input$start_analysis_TSa, {
  
  #    plot_df$species_freq <- plot_df$species_count/sum(plot_df$species_count)
  #    plot_df$seq_freq <- plot_df$seq_count/sum(plot_df$seq_count)
  #    plot_df$sample_freq <- plot_df$sample_count/sum(plot_df$sample_count)
  #    plot_df$`log2 fold change` <- log2(plot_df$seq_freq) - log2(plot_df$species_freq)
  #    plot_df$Labels <- paste0(plot_df$name, " (", round(plot_df$`log2 fold change\n(database spp - known spp)`, digits = 2), ")")
  
  plot_df$species_freq <- plot_df$species_count/sum(plot_df$species_count)
  plot_df$seq_freq <- plot_df$seq_count/sum(plot_df$seq_count)
  plot_df$sample_freq <- plot_df$sample_count/sum(plot_df$sample_count)
  plot_df$log2_fold_change <- round(log2(plot_df$seq_freq) - log2(plot_df$species_freq), digits = 2)
  plot_df$log2_fold_change_2 <- round(log2(plot_df$sample_freq) - log2(plot_df$species_freq), digits = 2)
  plot_df$color <- "gray"
  plot_df$color[plot_df$log2_fold_change > input$upper_bound] <- "orange"
  plot_df$color[plot_df$log2_fold_change < input$lower_bound] <- "blue"
  plot_df_regression <- plot_df[(plot_df$log2_fold_change > input$lower_bound & plot_df$log2_fold_change < input$upper_bound), ]
  
  #    plot_df$Labels <- paste0(plot_df$name, " (", round(plot_df$`log2 fold change`, digits = 2), ")")
  
  plot <- plot_df %>%
    
    ggplot(aes(x=species_freq, y=sample_freq)) +
    stat_cor(
      method = "spearman",
      label.x.npc = "left",
      label.y.npc = "top"
    ) +
    
#    geom_smooth(method=lm , color="gray", se=TRUE) + 
    
    geom_point(aes(name=name, log2_fold_change=log2_fold_change, log2_fold_change_2=log2_fold_change_2, sample_count=sample_count, size=species_count, color=color)) + 
    #                 text= paste("Name: ", Labels, "<br>", 
    #                             "Known species count: ", species_count, "<br>",
    #                             "Known species frequency: ", species_freq, "<br>",
    #                             "DB species frequency: ", seq_freq,
    #                             "log2 fold change: ", `log2 fold change`,
    #                             sep = ""))+
    
    scale_color_identity() +
    
    #      scale_color_gradient2(low = "#0096FF", mid = "#F7F7F7", high ="#F1A340", 
    #                           midpoint = 0, space = "Lab", guide = "colourbar") +
    
    scale_size_continuous(name = "Number of\nknown species",
                          breaks = c(1, 100, 1000, 10000, 25000, 50000),
                          limits = c(0, max(plot_df$species_count))) + 
    
    scale_x_log10(labels = scales::scientific) +
    scale_y_log10(labels = scales::scientific) +
    
    xlab("Known species") +
    ylab("Sampled species") +
    
    theme_light() +
    theme(
      axis.text=element_text(size=rel(1)),
      axis.title=element_text(size=rel(1))
    )
  
  plot <- plot + geom_smooth(data=plot_df_regression, aes(x=species_freq, y=seq_freq), method=lm , color="red", se=TRUE)
  
  ggplotly(plot, tooltip = c("name", "log2_fold_change", "log2_fold_change_2", "sample_count", "species_count"))
})

  output$get_scatterplot_DB_TS <- renderPlotly({ #renderPlotly?
    plot <- get_scatterplot_DB_TS()
    plot
  })


  output$get_scatterplot_DB <- renderPlotly({ #renderPlotly?
      plot <- get_scatterplot_DB()
      plot
  })

  output$get_scatterplot_TS <- renderPlotly({ #renderPlotly?
    plot <- get_scatterplot_TS()
    plot
  })
  
#  output$TS_data <- renderDataTable({
#    if (is.integer(output$count_data)) {
#    } else {
#      tmp <- get_TS_output()
#      tmp$spp_df
#    }
#  })
  
get_summary_count_db <- function(counts, root_taxon) {
  ancestry <- get_taxID_spp_counts(nodes = counts$nodes,
                                   start_from_species = FALSE,
                                   what = "ancestry",
                                   verbose = FALSE)
#  print(ancestry)
  
#  print(counts$countIDs)
  
  # Isolate only ancestries of sampled taxIDs
  idx <- which(sapply(ancestry,
                      function(x, ids){x[[1]] %in% ids},
                      ids = counts$countIDs))

#  print(idx)
    
  # Prepare nodes for joining
  nodes <- counts$nodes[, -2]
  names(nodes)[1] <- "taxID"
  
#  print(nodes)
  
  # Get sequence counts
  seqs <- data.table::data.table(taxID     = as.integer(names(counts$countIDs)),
                                 seq_count = as.integer(counts$countIDs))
  
#  print(seqs)
  
  # Get species counts
  spps <- data.table::as.data.table(counts$spp_df)

#  print(spps)
    
  # Join all info
  out <- seqs[spps, on = "taxID"]
  #out <- samples[out, on = "taxID"]
  out <- nodes[out, on = "taxID"]
#  print(out)
  idx <- which(sapply(ancestry,
                      function(x, tax){(tax %in% x) && (x[[1]] != tax)},
                      tax = root_taxon))
  
#  print(idx)
  
  ids <- sapply(ancestry[idx], function(x){x[[1]]})
  out <- out[out$taxID %in% ids, ]
  out
  #out$sample_count <- ifelse(is.na(out$sample_count), 0, out$sample_count)
}
  
})