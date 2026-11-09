library(shiny)
library(shinyFiles)
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


#if (require(bslib)) {
#  theme <- bslib::bs_theme(version = 4)
#} else {
#  theme <- NULL
#}


#holds the different parameters for user-defined tasks
parameter_tabs <- tabsetPanel(
  id = "params",
  type = "hidden",
  tabPanel("CreateDB",
             #Main title
             tags$h5("Create local NCBI Taxonomy count database\n"),
             tags$hr(),
             tags$h5("Parameters:\n"),
             tags$hr(),
           #Compute count database, or load a previously computed:
             tags$p('Path to user-defined tabular file linking NCBI TaxonIDs to sequence IDs:\n'),
             tags$p(),
             shinyFilesButton("user_database_file_path", "File select", "Please select a file", multiple = TRUE, viewtype = "detail"),
             tags$hr(),
             tags$p("Path to folder containing NCBI Taxonomy dump (.dmp) files:\n"),
             tags$p(),
             shinyDirButton("NCBI_Taxonomy_directory", "Folder select", "Please select a folder"),
             tags$hr(),
             tags$p("Path to folder to save database output files:\n"),
             tags$p(),
             shinyDirButton("TS_DB_database_path", "Folder select", "Please select a folder"),
             tags$hr(),
             actionButton("start_analysis_DB_build", "Build TaxonSampling Database"),
  ),
  tabPanel("AnalyzeDB", 
             #Main title
           tags$h5("Analyze previously extracted database\n"),
           tags$hr(),
           tags$h5("Parameters:\n"),
           tags$hr(),
           tags$p("Path to previously built TaxonSampling database:\n"),
             tags$hr(),
           shinyDirButton("TS_DB_database_path_previous_run", "Folder select", "Please select a folder"),
           tags$hr(),
           selectInput("taxonomic_level",
                         tags$p("Taxonomic level to survey database:"), 
                         choices = list(
                                        #                                     "Subgenus" = "subgenus",
                                        #                                     "No rank" = "no rank",
                                        #                                     "Species group" = "species group",
                                        #                                     "Subfamily" = "subfamily",
                                        "Phylum" = "phylum",
                                        "Class" = "class",
                                        "Order" = "order",
                                        "Family" = "family",
                                        "Genus" = "genus",
                                        "Species" = "species"
                                        #                                     "Clade" = "clade",
                                        #                                     "Suborder" = "suborder",
                                        #                                     "Tribe" = "tribe",
                                        #                                     "Superfamily" = "superfamily",
                                        #                                     "Parvorder" = "parvorder",
                                        #                                     "Infraclass" = "infraclass",
                                        #                                     "Infraorder" = "infraorder",
                                        #                                     "Subclass" = "subclass",
                                        #                                     "Superorder" = "superorder",
                                        #                                     "Superclass" = "superclass",
                                        #                                     "Subcohort" = "subcohort",
                                        #                                     "Cohort" = "cohort",
#                                        "Subspecies" = "subspecies"
                                        
                         ), selected = "class"),
             tags$hr(),
             numericInput("root_taxon", 
                          tags$p("NCBI TaxID for root taxon:"), 
                          value = 2759),   
             tags$hr(),
             numericInput("upper_bound", 
                          tags$p("Upper bound for regression:"), 
                          value = 0.5),   
             tags$hr(),
             numericInput("lower_bound", 
                          tags$p("Lower bound for regression:"), 
                          value = -0.5),   
             tags$hr(),
             actionButton("start_analysis_DB_plot", "Analyze DB"),
  ),
  tabPanel("CreateTS",
           #Main title
           tags$h5("Sample taxa from database using TaxonSampling\n"),
           tags$hr(),
           tags$h5("Parameters:\n"),
           tags$hr(),
           tags$p("Path to a previously built TaxonSampling database:\n"),
           shinyDirButton("TS_DB_database_input_path", "Folder select", "Please select a folder"),
           tags$hr(),
           tags$p("Path to directory to write TS output:\n"),
           shinyDirButton("TS_outptut_path", "Folder select", "Please select a folder"),
           tags$hr(),
           numericInput("m",
             tags$p("Number of terminal nodes (species) to sample:"), 
             value = 500),
           tags$hr(),
           numericInput("root_taxon", 
                       tags$p("NCBI TaxID for root taxon:"), 
                       value = 2759),
           tags$hr(),
           selectInput("sampling",
            tags$p("Sampling frequency:"), 
            choices = list("Agnostic" = "agnostic", "Known species" = "known_species"
                           ), selected = "known_species"),          
           tags$hr(),
           selectInput("method",
            tags$p("Sampling method:"), 
            choices = list("Balanced" = "balanced", "Diversity" = "diversity"
            ), selected = "diversity"),
           tags$hr(),
           selectInput("randomize",
            tags$p("Randomize?"), 
            choices = list("Yes" = "yes", "No" = "no", "After first round" = "after_first_round"
            ), selected = "after_first_round"),
          actionButton("start_analysis_TS", "Run TaxonSampling"),
          tags$hr(),
  ),

  tabPanel("AnalyzeTS",
           tags$h5("Compare DB and TaxonSampling data\n"),
           tags$hr(),
           tags$h5("Parameters:\n"),
           tags$hr(),
           tags$p("Path to the directory of a previously extracted database:\n"),
           shinyDirButton("TS_DB_input_path", "Folder select", "Please select a folder"),
           tags$hr(),
           tags$p("Path to the directory of a previous TaxonSampling run:\n"),
           shinyDirButton("TS_run_input_path", "Folder select", "Please select a folder"),
           tags$hr(),
           selectInput("taxonomic_level_TS",
                       tags$p("Taxonomic level to survey DB/TS:"), 
                       choices = list(
                         #                                     "Subgenus" = "subgenus",
                         #                                     "No rank" = "no rank",
                         #                                     "Species group" = "species group",
                         #                                     "Subfamily" = "subfamily",
                         "Phylum" = "phylum",
                         "Class" = "class",
                         "Order" = "order",
                         "Family" = "family",
                         "Genus" = "genus",
                         "Species" = "species"
                         #                                     "Clade" = "clade",
                         #                                     "Suborder" = "suborder",
                         #                                     "Tribe" = "tribe",
                         #                                     "Superfamily" = "superfamily",
                         #                                     "Parvorder" = "parvorder",
                         #                                     "Infraclass" = "infraclass",
                         #                                     "Infraorder" = "infraorder",
                         #                                     "Subclass" = "subclass",
                         #                                     "Superorder" = "superorder",
                         #                                     "Superclass" = "superclass",
                         #                                     "Subcohort" = "subcohort",
                         #                                     "Cohort" = "cohort",
                         #                                        "Subspecies" = "subspecies"
                         
                       ), selected = "class"),
           tags$hr(),
           numericInput("upper_bound_TS", 
                        tags$p("Upper bound for regression:"), 
                        value = 0.5),   
           tags$hr(),
           numericInput("lower_bound_TS", 
                        tags$p("Lower bound for regression:"), 
                        value = -0.5),   
           
           tags$hr(),
           actionButton("start_analysis_TS_plot", "Analyze TaxonSampling result"),
  )
)

fluidPage(
  theme = theme, 
  img(src = "TS_logo.png",
      height="8%",
      width="8%"
      ),

#  headerPanel(
#    "TaxonSampling"
#  ),
  sidebarLayout(
    sidebarPanel(
      selectInput("task", tags$h4("Task selection:"),
                  choices = c("CreateDB", "AnalyzeDB", "CreateTS", "AnalyzeTS")
      ),
      parameter_tabs,
    ),
#  ),

    mainPanel(
      tabsetPanel(id="outtabs",
                  type="tabs",
        tabPanel(title="Create DB",
                 value="CreateDB",
          tags$h4("Input file paths"),
          tags$hr(),
          tags$p(HTML("Path to local TaxID2SeqID file:\n")),
          verbatimTextOutput("TaxID2SeqID_file_path"),
#          tags$hr(),
          tags$p(HTML("Path to local NCBI Taxonomy database:\n")),
          verbatimTextOutput("NCBI_Tax_path"),
#          tags$hr(),
          tags$p(HTML("Path to write output R database object:\n")),
          verbatimTextOutput("saveDir_path"),
          tags$hr(),
          tags$h4("Execution log"),
          tags$hr(),
          verbatimTextOutput("create_countsDB"),
        ),
        tabPanel(title="DB summary",
                 value="AnalyzeDB",
        tags$h4("Database analysis results"),
        tags$hr(),
        tags$h5("Scatterplot of known species (x) vs species in databse (y)"),
        plotlyOutput("get_scatterplot_DB"),
        tags$hr(),
        tags$h5("Taxonomic diversity and biases found in database"),
        DT::dataTableOutput("DB_summary"),
        tags$hr(),
        ),
        tabPanel(title="Create TS",
                 value="CreateTS",
          tags$h4("Input file paths"),
          tags$hr(),
          tags$p(HTML("Path to local TaxonSampling DB:\n")),
          verbatimTextOutput("inputDir_TS_path"),
          tags$p(HTML("Path to write TaxonSampling output:\n")),
          verbatimTextOutput("saveDir_TS_path"),
          tags$hr(),
          tags$h4("Execution log"),
          tags$hr(),
          verbatimTextOutput("create_TSDB"),
#          tableOutput("DB_summary"),
#          tags$p(HTML('Database profile:')),
#          tags$p("Database description:"),
#          plotlyOutput("get_scatterplot_DB_TS"),
#          tags$p("TaxonSampling result:"),
#          tags$p(HTML('TaxonSampling profile:')),
#          plotlyOutput("get_scatterplot_TS")
        ),

        tabPanel(title="TS summary",
                 value="AnalyzeTS",
          tags$h4("Input file paths"),
          tags$hr(),
          tags$p(HTML("Path to directory with DB output:")),
          verbatimTextOutput("TS_DB_input_path"),
          tags$p(HTML("Path directory with TS output:")),
          verbatimTextOutput("TS_run_input_path"),
#          verbatimTextOutput("create_TSDB"),
          tags$p(HTML('Database profile:')),
          plotlyOutput("get_scatterplot_DB_TS"),
          tags$p(HTML('TaxonSampling profile:')),
          plotlyOutput("get_scatterplot_TS"),
          tags$p(HTML("Taxonomic representativity in TS sampling:")),
          DT::dataTableOutput("TS_summary"),
          tags$p(),
          tags$hr(),
        )
      )
    )
  )
)
