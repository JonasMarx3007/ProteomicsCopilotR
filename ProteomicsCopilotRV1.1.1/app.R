#V0.1.1
library(shiny)
library(readr)
library(readxl)
library(dplyr)
library(reshape2)
library(ggplot2)
library(stringr)
library(plotly)
library(tidyverse)
library(gprofiler2)
library(corrplot)
library(pheatmap)
library(reticulate)
library(xtable)
library(ggeasy)
library(tools)
library(zip)
library(DT)
library(rsvg)
library(png)
library(protr)
library(r3dmol)
library(UniprotR)
library(protti)
library(KSEAapp)
library(shinycssloaders)
library(svglite)
library(rstudioapi)
library(seqinr)
library(openxlsx)
options(warn = -1)
counter_file <- "www/ressources/logcounter.txt"
workmode = "funn"

source("functions.R")

#### UI functions ####
ui <- fluidPage(
  title="Proteomics Copilot",
  
  tags$head(
    tags$link(rel = "shortcut icon", href = "icon/receptor_icon.ico")
  ),
  
  tags$head(
    tags$style(HTML("
      /* Tab background */
      .nav-tabs {
        background-color: #337ab7;
        border-bottom: 1px solid #2a5d8f;
      }
      /* Inactive tabs */
      .nav-tabs > li > a {
        color: white;
        background-color: #337ab7;
        border: 1px solid #2a5d8f;
        border-bottom-color: transparent;
      }
      /* Hover state */
      .nav-tabs > li > a:hover {
        background-color: #2a5d8f;
      }
      /* Active tab */
      .nav-tabs > li.active > a,
      .nav-tabs > li.active > a:focus {
        color: white;
        background-color: #f39c12;
        border: 1px solid #d78e10;
        border-bottom-color: #f39c12;
      }
    "))
  ),
  
  titlePanel(div(
    "Proteomics Copilot",
    style = "
        background-color: #337ab7;
        color: white;
        padding: 10px;
        border-radius: 4px;
      "
  )),
  fluidRow(
    column(width = 12,
           tabsetPanel(
#### SuperTab01 - Data ####
tabPanel("Data",
          tabsetPanel(
#### Tab01 - Data Upload ####
              tabPanel(
                "Data Upload",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    h4("Upload Files"),
                    fileInput("file", "Protein Level Data", 
                              accept = c(".csv", ".tsv", ".txt", ".xlsx")),
                    fileInput("file3", "Phospho Data", 
                              accept = c(".csv", ".tsv", ".txt", ".xlsx")),
                    fileInput("file2", "Full Report", 
                              accept = c(".csv", ".tsv", ".txt", ".xlsx")),
                    hr(),
                    h4("Collapse Options"),
                    radioButtons("collapse_option", "Select Collapse Option:",
                                 choices = c("Collapsed" = "collapsed", 
                                             "Not Collapsed" = "not_collapsed"),
                                 selected = "collapsed"),
                    numericInput("collapse_cutoff", "Collapse Cutoff:", value = 0)
                  ),
                  mainPanel(
                    h4("Protein Level Data Preview"),
                    DT::dataTableOutput("table"),
                    hr(),
                    h4("Phospho Data Preview"),
                    DT::dataTableOutput("table3"),
                    hr(),
                    h4("Full Report Preview"),
                    DT::dataTableOutput("table2")
                  )
                )
              ),
#### Tab02 - Meta annotation #####
            tabPanel(
              "Data Annotation",
              sidebarLayout(
                sidebarPanel(
                  width = 3,
                  
                  h4("Normal Data Annotation"),
                  fileInput("upload_meta", "Upload Metadata (Protein Group)", accept = c(".csv", ".xlsx", ".txt", ".tsv")),
                  h5("Is the data log2 transformed?"),
                  fluidRow(
                    column(6, actionButton("log2_yes", "Yes")),
                    column(6, actionButton("log2_no", "No"))
                  ),
                  hr(),
                  numericInput("filter_num", "Filter: At least", value = 3, min = 1),
                  selectInput("filterop1", "Value(s)", choices = c("per group", "in at least one group")),
                  actionButton("apply_filter", "Apply Filter"),
                  
                  hr(),
                  
                  h4("Phospho Data Annotation"),
                  fileInput("upload_meta2", "Upload Metadata (Phospho)", accept = c(".csv", ".xlsx", ".txt", ".tsv")),
                  h5("Is the data log2 transformed?"),
                  fluidRow(
                    column(6, actionButton("log2_yes2", "Yes")),
                    column(6, actionButton("log2_no2", "No"))
                  ),
                  hr(),
                  numericInput("filter_num2", "Filter: At least", value = 3, min = 1),
                  selectInput("filterop2", "Value(s)", choices = c("per group", "in at least one group")),
                  actionButton("apply_filter2", "Apply Filter"),
                  
                  hr(),
                  
                  h4("Color Scheme"),
                  selectInput(
                    inputId = "color_palette",
                    label = "Choose a Color Palette:",
                    choices = c("Default", "Mario Document Input", "Default16", "Warm/Cold", "Black/Grey", "Yue7"),
                    selected = "Default"
                  ),
                  actionButton("reloadButton", "Reload All Plots")
                ),
                
                mainPanel(
                  h4("Normal Data Condition Setup"),
                  numericInput("num_conditions", "Number of Conditions:", value = 1, min = 1),
                  uiOutput("condition_inputs"),
                  hr(),
                  
                  h4("Annotated Normal Data"),
                  DT::dataTableOutput("displayed_data"),
                  hr(),
                  
                  h4("Phospho Data Condition Setup"),
                  numericInput("num_conditions2", "Number of Conditions:", value = 1, min = 1),
                  uiOutput("condition_inputs2"),
                  hr(),
                  
                  h4("Annotated Phospho Data"),
                  DT::dataTableOutput("displayed_data2")
                )
              )
            ),
#### Tabn1 - Data imputation ####
            tabPanel(
              "Impute Data",
              sidebarLayout(
                sidebarPanel(
                  width = 3,
                  actionButton("ImputeEVE", "Impute Data Values"),
                  hr(),
                  h4("Imputation Settings"),
                  numericInput("qn1", "q-Value:", value = 0.01),
                  numericInput("adj_stdn1", "Adjust Standard Deviation:", value = 1),
                  numericInput("seedn1", "Random Seed:", value = 1337),
                  selectInput("leveln1", 
                              "Level:", 
                              choices = c("Protein", "Phosphosite")),
                  hr(),
                  downloadButton("imputed_data_down", "Download Imputed Data")
                ),
                mainPanel(
                  h4("Imputation Plots"),
                  customSpinnerWrapper(plotOutput("Impute1"), workmode = workmode),
                  customSpinnerWrapper(plotOutput("Impute2"), workmode = workmode),
                  customSpinnerWrapper(plotOutput("Impute3"), workmode = workmode),
                  h4("Imputed Data Table"),
                  DT::dataTableOutput("imputed_data_tab")
                )
              )
            ),
#### Tab4.5 - Data Distribution ####
            tabPanel(
              "Distribution",
              sidebarLayout(
                sidebarPanel(
                  width = 3,
                  selectInput("level4.5", "Level:",
                              choices = c("Protein", "Phosphosite"))
                ),
                mainPanel(
                  customSpinnerWrapper(
                    plotOutput("qqnorm", height = "600px"),
                    workmode = workmode
                  ),
                  tags$img(
                    src = "material/qqnorm.jpg", 
                    width = "100%", 
                    style = "margin-bottom: 20px;"
                  ),
                  tags$img(
                    src = "material/qqnorm_txt.jpg", 
                    width = "100%", 
                    style = "margin-bottom: 20px;"
                  )
                )
              )
            ),

#### Tab3.5 - Data Verification ####
          tabPanel(
            "Verification",
            sidebarLayout(
              sidebarPanel(
                width = 3,
                selectInput("level3.5", "Level:",
                            choices = c("Protein", "Phosphosite")),
              ),
              mainPanel(
                customSpinnerWrapper(
                  plotOutput("FirstDigitPlot", height = "600px"),
                  workmode = workmode
                ),
                customSpinnerWrapper(
                  plotOutput("DataStructurePlot", height = "600px"),
                  workmode = workmode
                )
              )
            )
          ),
)),
#### SuperTab02 - QC Pipeline #####
tabPanel("QC Pipeline",
          tabsetPanel(
#### Tab03 - Coverage Plot ####
              tabPanel(
                "Coverage Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("toggle_id3", "Toggle IDs"),
                    actionButton("toggle_header3", "Toggle Header"),
                    actionButton("toggle_legend3", "Toggle Legend"),
                    
                    hr(),
                    selectInput("level3", "Level:",
                                choices = c("Protein", "Peptide", "Phosphosite")),
                    hr(),
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth3", "Width (cm):", value = 20),
                    numericInput("plotHeight3", "Height (cm):", value = 10),
                    numericInput("plotDPI3", "DPI:", value = 300),
                    selectInput("plotFormat3", "File Format:",
                                choices = c("PNG" = "png", "JPG" = "jpg", "SVG" = "svg", "PDF" = "pdf"),
                                selected = "png"),
                    downloadButton("downloadCoveragePlot", "Download Plot"),
                    hr(),
                    
                    selectInput("textPosition3", "Position:",
                                choices = c("Above" = "up", "Below" = "down"),
                                selected = "up"),
                    
                    textAreaInput("text3", "Annotation Text:", value = "",
                                  width = '100%', height = '100px'),
                    
                    actionButton("addText3", "Add"),
                    actionButton("deleteText3", "Delete")
                  ),
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("coveragePlot", height = "600px"),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab04 - Missing Values Plot ####
              tabPanel(
                "Missing Value Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    actionButton("toggle_header4", "Toggle Header"),
                    hr(),
                    
                    selectInput("level4", "Level:",
                                choices = c("Protein", "Peptide", "Phosphosite", "Precursor")),
                    
                    numericInput("missValBin4", "Bin missing values (optional):", value = 0, min = 0, step = 1),
                    hr(),
                    
                    numericInput("text_size4",  "Text Size:",  value = 3.88),
                    actionButton("toggle_text4", "Toggle Text"),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth4", "Width (cm):", value = 20),
                    numericInput("plotHeight4", "Height (cm):", value = 10),
                    numericInput("plotDPI4", "DPI:", value = 300),
                    
                    selectInput("missValPlotFormat", "Download Format:",
                                choices = c("PNG" = "png", "JPG" = "jpg", "SVG" = "svg", "PDF" = "pdf"),
                                selected = "png"),
                    
                    downloadButton("downloadMissValPlot", "Download Plot"),
                    
                    hr(),
                    
                    selectInput("textPosition4", "Position:",
                                choices = c("Above" = "up", "Below" = "down"),
                                selected = "up"),
                    
                    textAreaInput("text4", "Annotation Text:", value = "",
                                  width = '100%', height = '100px'),
                    
                    actionButton("addText4", "Add"),
                    actionButton("deleteText4", "Delete")
                  ),
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("MissValPlot", height = "600px"),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab05 - Histogramm ####
              tabPanel(
                "Histogram Intensity",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    actionButton("toggle_header5", "Toggle Header"),  
                    actionButton("toggle_legend5", "Toggle Legend"),   
                    hr(),
                    selectInput("level5", "Level:",
                                choices = c("Protein", "Phosphosite")),
                    hr(),
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth5", "Width (cm):", value = 20),
                    numericInput("plotHeight5", "Height (cm):", value = 10),
                    numericInput("plotDPI5", "DPI:", value = 300),
                    selectInput("plotFormat5", "File Format:",
                                choices = c("PNG" = "png", "JPG" = "jpg", "SVG" = "svg", "PDF" = "pdf"),
                                selected = "png"),
                    downloadButton("downloadHistIntPlot", "Download Plot"),
                    hr(),
                    
                    selectInput("textPosition5", "Position:",
                                choices = c("Above" = "up", "Below" = "down"),
                                selected = "up"),
                    
                    textAreaInput("text5", "Annotation Text:", value = "",
                                  width = '100%', height = '100px'),
                    
                    actionButton("addText5", "Add"),
                    actionButton("deleteText5", "Delete")
                  ),
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("HistIntPlot", height = "600px"),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab06 - Boxplot ####
              tabPanel(
                "Boxplot Intensity",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("toggle_outliersBX", "Toggle Outliers"),
                    actionButton("toggle_meanBX", "Mean/Single"),
                    actionButton("toggle_id6", "Toggle ID"),
                    actionButton("toggle_header6", "Toggle Header"),  
                    actionButton("toggle_legend6", "Toggle Legend"),   
                    hr(),
                    
                    selectInput("level6", "Level:", choices = c("Protein", "Phosphosite")),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth6", "Width (cm):", value = 20),
                    numericInput("plotHeight6", "Height (cm):", value = 10),
                    numericInput("plotDPI6", "DPI:", value = 300),
                    selectInput("plotFormat6", "File Format:",
                                choices = c("PNG" = "png", "JPG" = "jpg", "SVG" = "svg", "PDF" = "pdf"),
                                selected = "png"),
                    downloadButton("downloadBoxIntPlot", "Download Plot"),
                    hr(),
                    
                    selectInput("textPosition6", "Position:",
                                choices = c("Above" = "up", "Below" = "down"),
                                selected = "up"),
                    
                    textAreaInput("text6", "Annotation Text:", value = "",
                                  width = '100%', height = '100px'),
                    
                    actionButton("addText6", "Add"),
                    actionButton("deleteText6", "Delete")
                  ),
                  
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("BoxIntPlot", height = "600px"),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab07 - COV Plot ####
              tabPanel(
                "Coefficient of Variation Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    actionButton("toggle_outliersCOV", "Toggle Outliers"),
                    actionButton("toggle_header7", "Toggle Header"),
                    actionButton("toggle_legend7", "Toggle Legend"),
                    hr(),
                    selectInput("level7", "Level:", choices = c("Protein", "Phosphosite")),
                    hr(),
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth7", "Width (cm):", value = 20),
                    numericInput("plotHeight7", "Height (cm):", value = 10),
                    numericInput("plotDPI7", "DPI:", value = 300),
                    selectInput("plotFormat7", "Format:", choices = c("PNG" = "png", "JPG" = "jpg", "SVG" = "svg", "PDF" = "pdf")),
                    downloadButton("downloadCovPlot", "Download Plot"),
                    hr(),
                    
                    selectInput("textPosition7", "Position:",
                                choices = c("Above" = "up", "Below" = "down"),
                                selected = "up"),
                    
                    textAreaInput("text7", "Annotation Text:", value = "", width = '100%', height = '100px'),
                    
                    actionButton("addText7", "Add"),
                    actionButton("deleteText7", "Delete")
                  ),
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("CovPlot", height = "600px"),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab08 - PCA #### 
              tabPanel(
                "PCA Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    actionButton("toggle_header8", "Toggle Header"),
                    actionButton("toggle_legend8", "Toggle Legend"),
                    hr(),
                    selectInput("level8", "Level:", choices = c("Protein", "Phosphosite")),
                    selectInput("style8", "Select Plot Type:", choices = c("PCA", "tSNE", "UMAP")),
                    hr(),
                    selectInput("type8", "Type:", choices = c("Non-Interactive", "Interactive"), selected = "Non-Interactive"),
                    sliderInput("dotSize8", "Dot Size:", min = 1, max = 20, value = 3),
                    hr(),
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth8", "Width (cm):", value = 20),
                    numericInput("plotHeight8", "Height (cm):", value = 10),
                    numericInput("plotDPI8", "DPI:", value = 300),
                    selectInput("plotFormat8", "File Format:", choices = c("PNG" = "png", "JPG" = "jpg", "SVG" = "svg", "PDF" = "pdf")),
                    downloadButton("downloadPCAPlot", "Download Plot"),
                    hr(),
                    selectInput("textPosition8", "Position:",
                                choices = c("Above" = "up", "Below" = "down"),
                                selected = "up"),
                    textAreaInput("text8", "Annotation Text:", value = "", width = '100%', height = '100px'),
                    actionButton("addText8", "Add"),
                    actionButton("deleteText8", "Delete")
                  ),
                  mainPanel(
                    customSpinnerWrapper(
                      tagList(
                        conditionalPanel(
                          condition = "input.type8 == 'Non-Interactive'",
                          plotOutput("PCAPlot", height = "600px")
                        ),
                        conditionalPanel(
                          condition = "input.type8 == 'Interactive'",
                          plotlyOutput("PCAPlotInteractive", height = "600px")
                        )
                      ),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab09 - Abundance Plot ####
              tabPanel(
                "Abundance Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    selectInput("level9", "Level:", choices = c("Protein", "Phosphosite")),
                    hr(),
                    h4("Plot Size & Download"),
                    numericInput("plotWidth9", "Width (cm):", value = 20),
                    numericInput("plotHeight9", "Height (cm):", value = 10),
                    numericInput("plotDPI9", "DPI:", value = 300),
                    downloadButton("downloadAbPlot", "Download Plot"),
                    hr(),
                    textAreaInput("text9", "Enter text:", value = "", width = '100%', height = '150px'),
                    selectInput("textPosition9", "Position:", choices = c("Above" = "up", "Below" = "down"), selected = "up"),
                    actionButton("addText9", "Add"),
                    actionButton("deleteText9", "Delete")
                  ),
                  mainPanel(
                    fluidRow(
                      column(
                        width = 6,
                        selectInput("condition9", "Choose Condition:", choices = c("All Conditions"), selected = "All Conditions")
                      ),
                      column(
                        width = 6,
                        selectizeInput("protein9", "Select Proteins:", choices = NULL, multiple = TRUE, width = '100%')
                      )
                    ),
                    hr(),
                    customSpinnerWrapper(
                      plotlyOutput("abundancePlot", height = "600px"),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab12 - Correlation Plot ####
              tabPanel(
                "Correlation Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("Change12", "Change Display"),
                    actionButton("toggle_id12", "Toggle ID"),
                    hr(),
                    
                    selectInput("level12", "Level:", choices = c("Protein", "Phosphosite")),
                    hr(),
                    
                    selectInput("textPosition12", "Position:", choices = c("Above" = "up", "Below" = "down"), selected = "up"),
                    
                    textAreaInput("text12", "Annotation Text:", value = "", width = '100%', height = '100px'),
                    
                    actionButton("addText12", "Add"),
                    actionButton("deleteText12", "Delete")
                  ),
                  
                  mainPanel(
                    h4("Correlation Plot (Pearson)"),
                    customSpinnerWrapper(
                      plotOutput("CorrPlot", height = "600px"),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab13 - Heatmap ####
              tabPanel(
                "Heatmap",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("toggle_id13", "Toggle ID"),
                    actionButton("show13", "Show Missing Values"),
                    actionButton("fix13", "Fix Plot"),
                    hr(),
                    
                    selectInput("level13", "Level:", choices = c("Protein", "Phosphosite")),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth13", "Width (cm):", value = 20),
                    numericInput("plotHeight13", "Height (cm):", value = 10),
                    numericInput("plotDPI13", "DPI:", value = 300),
                    hr(),
                    
                    downloadButton("downloadheatPlot", "Download Plot"),
                    hr(),
                    
                    selectInput("textPosition13", "Position:", choices = c("Above" = "up", "Below" = "down"), selected = "up"),
                    textAreaInput("text13", "Annotation Text:", value = "", width = "100%", height = "100px"),
                    actionButton("addText13", "Add"),
                    actionButton("deleteText13", "Delete")
                  ),
                  
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("heatPlot", height = "800px"),
                      workmode = workmode
                    )
                  )
                )
              ),
)),
#### SuberTab03 - Statistical Analysis ####
tabPanel("Statistical Analysis",
          tabsetPanel(
#### Tab10 - Volcano Plot ####
              tabPanel(
                "Volcano Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    h4("Thresholds"),
                    numericInput("in_pval10", "P‑value threshold:", value = 0.05, step = 0.01),
                    numericInput("in_log2fc10", "log₂ FC threshold:", value = 1, step = 0.1),
                    checkboxInput("uncorrected10", "Use uncorrected p‑values", value = FALSE),
                    hr(),
                    
                    selectInput("paired10", "Test Type:", choices = c("Unpaired", "Paired")),
                    selectInput("level10", "Level:", choices = c("Protein", "Phosphosite")),
                    hr(),
                    
                    downloadButton("download10", "Download Volcano Data"),
                    downloadButton("download10a", "Download All Volcano Data"),
                    hr(),
                    
                    selectInput("textPosition10", "Position:", choices = c("Above" = "up", "Below" = "down"), selected = "up"),
                    textAreaInput("text10", "Annotation Text:", value = "", width = '100%', height = '100px'),
                    actionButton("addText10", "Add"),
                    actionButton("deleteText10", "Delete"),
                    hr(),
                    
                    actionButton("addVolc", "Add Plot to Report", 
                                 icon = icon("file-medical"), 
                                 style = "width:100%; margin-top: 20px;"),
                    
                    actionButton("removeVolc", "Remove Plot from Report", 
                                 icon = icon("trash"), 
                                 style = "width:100%; margin-top: 10px; background-color: #dc3545; color: white;")
                  ),
                  
                  mainPanel(
                    fluidRow(
                      column(6, selectInput("condition1_10", "Condition 1:", choices = c())),
                      column(6, selectInput("condition2_10", "Condition 2:", choices = c()))
                    ),
                    hr(),
                    
                    customSpinnerWrapper(
                      plotlyOutput("VolcPlot", height = "500px"),
                      workmode = workmode
                    ),
                    
                    hr(),
                    DT::dataTableOutput("table10")
                  )
                )
              ),
#### Tab11 - GSEA ####
              tabPanel(
                "GSEA",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    numericInput(
                      "top_n11",
                      "Number of top gene sets to show:",
                      value = 10,
                      min   = 1,
                      step  = 1
                    ),
                    hr(),
                    
                    h4("Term Size Filter"),
                    numericInput(
                      "filter11min",
                      "Min term size:",
                      value = 20
                    ),
                    numericInput(
                      "filter11max",
                      "Max term size:",
                      value = 300
                    ),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput(
                      "plotWidth11",
                      "Width (cm):",
                      value = 20
                    ),
                    numericInput(
                      "plotHeight11",
                      "Height (cm):",
                      value = 10
                    ),
                    numericInput(
                      "plotDPI11",
                      "DPI:",
                      value = 300
                    ),
                    hr(),
                    
                    downloadButton("downloadUpPlot",   "Download Plot (UP)"),
                    downloadButton("downloadDownPlot", "Download Plot (DOWN)")
                  ),
                  
                  mainPanel(
                    h3("Upregulated Gene Sets"),
                    customSpinnerWrapper(
                      plotOutput("UpregEnrichmentPlot", height = "600px"),
                      workmode=workmode
                    ),
                    hr(),
                    
                    h3("Downregulated Gene Sets"),
                    customSpinnerWrapper(
                      plotOutput("DownregEnrichmentPlot", height = "600px"),
                      workmode=workmode
                    ),
                    hr(),
                    
                    fluidRow(
                      column(
                        6,
                        h4("Upregulated Gene List"),
                        DT::dataTableOutput("table11up")
                      ),
                      column(
                        6,
                        h4("Downregulated Gene List"),
                        DT::dataTableOutput("table11down")
                      )
                    )
                  )
                )
              ),
### Tab10.5 Simulation ####
              tabPanel(
                "Simulation",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    h4("Thresholds"),
                    numericInput("in_pval10.5",   "P‑value threshold:",    value = 0.05, step = 0.01),
                    numericInput("in_log2fc10.5", "log₂ FC threshold:",    value = 1,    step = 0.1),
                    hr(),
                    
                    selectInput("level10.5",  "Level:", choices = c("Protein", "Phosphosite")),
                    hr(),
                    
                    sliderInput("mod_var10.5", "Variance multiplier:", min = 0.25, max = 4, value = 1, step = 0.05),
                    sliderInput("mod_n10.5",   "Sample size override (n):", min = 1, max = 20, value = 1, step = 1)
                  ),
                  
                  mainPanel(
                    fluidRow(
                      column(6, selectInput("condition1_10.5", "Condition 1:", choices = c())),
                      column(6, selectInput("condition2_10.5", "Condition 2:", choices = c()))
                    ),
                    hr(),
                    
                    customSpinnerWrapper(
                      plotlyOutput("VolcPlotSim", height = "500px"),
                      workmode = workmode
                    )
                  )
                )
              ),
)),
#### SuperTab04 - Peptide Level #####
tabPanel("Peptide Level",
         tabsetPanel(
#### Tab14 - RT Plot ####
              tabPanel(
                "RT Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("line14", "Add Line"),
                    actionButton("toggle_header14", "Toggle Header"),
                    hr(),
                    
                    selectInput("style14", "Select Plot Type:",
                                choices = c("Scatter Plot", "Hexbin Plot", "Density Plot")),
                    numericInput("bins14", "Bins:", value = 1000),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth14", "Width (cm):", value = 20),
                    numericInput("plotHeight14", "Height (cm):", value = 10),
                    numericInput("plotDPI14", "DPI:", value = 300),
                    selectInput("plotFormat14", "File Format:",
                                choices = c("PNG" = "png", "JPG" = "jpg", "SVG" = "svg", "PDF" = "pdf"),
                                selected = "png"),
                    downloadButton("downloadRTPlot", "Download Plot"),
                    hr(),
                    
                    selectInput("textPosition14", "Position:", choices = c("Above" = "up", "Below" = "down"), selected = "up"),
                    textAreaInput("text14", "Annotation Text:", value = "", width = '100%', height = '100px'),
                    actionButton("addText14", "Add"),
                    actionButton("deleteText14", "Delete")
                  ),
                  
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("RTPlot", height = "600px"),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab15 - Modification Plot ####
              tabPanel(
                "Modification Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("toggle_id15", "Toggle ID"),
                    actionButton("toggle_header15", "Toggle Header"),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth15", "Width (cm):", value = 20),
                    numericInput("plotHeight15", "Height (cm):", value = 10),
                    numericInput("plotDPI15", "DPI:", value = 300),
                    selectInput("plotFormat15", "File Format:",
                                choices = c("PNG" = "png", "JPG" = "jpg", "SVG" = "svg", "PDF" = "pdf"),
                                selected = "png"),
                    downloadButton("downloadModPlot", "Download Plot"),
                    hr(),
                    
                    selectInput("textPosition15", "Position:", choices = c("Above" = "up", "Below" = "down")),
                    textAreaInput("text15", "Enter text:", value = "", width = "100%", height = "200px"),
                    actionButton("addText15", "Add"),
                    actionButton("deleteText15", "Delete")
                  ),
                  
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("ModPlot", height = "600px"),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab16 - Missed Cleavage Plot ####
              tabPanel(
                "Missed Cleavage Plot",
                sidebarLayout(
                  sidebarPanel(
                    width = 3,
                    
                    actionButton("toggle_id16", "Toggle ID"),
                    actionButton("toggle_header16", "Toggle Header"),
                    hr(),
                    
                    numericInput("text_size16", "Text Size:", value = 3.88),
                    actionButton("toggle_text16", "Toggle Text"),
                    hr(),
                    
                    h4("Plot Size & Resolution"),
                    numericInput("plotWidth16", "Width (cm):", value = 20),
                    numericInput("plotHeight16", "Height (cm):", value = 10),
                    numericInput("plotDPI16", "DPI:", value = 300),
                    selectInput("plotFormat16", "File Format:",
                                choices = c("PNG" = "png", "JPG" = "jpg", "SVG" = "svg", "PDF" = "pdf"),
                                selected = "png"),
                    downloadButton("downloadMCPlot", "Download Plot"),
                    hr(),
                    
                    selectInput("textPosition16", "Position:", choices = c("Above" = "up", "Below" = "down")),
                    textAreaInput("text16", "Enter text:", value = "", width = "100%", height = "200px"),
                    actionButton("addText16", "Add"),
                    actionButton("deleteText16", "Delete")
                  ),
                  
                  mainPanel(
                    customSpinnerWrapper(
                      plotOutput("MCPlot", height = "600px"),
                      workmode = workmode
                    )
                  )
                )
              ),
#### Tab22 - Protein centric ####
              tabPanel("Protein Centric",
                    actionButton("db_fasta22", "Load Database"),
                    verbatimTextOutput("db_status"),         
                    selectizeInput("protein20", "Select Protein:", choices = NULL),
                    numericInput("chunk20", "Chunk size:", value = 150),
                    verbatimTextOutput("seq_allign"),
                    r3dmolOutput("plot_3d22")
              ),
)),
#### SuberTab05 - Single Protein ####
tabPanel("Single Protein",
         tabsetPanel(
#### Tab17 - Boxplot, single protein ####
             tabPanel("Protein Box",
                      selectizeInput("protein17", "Select Protein:", choices = NULL),
                      selectizeInput("conditions17", "Select Conditions:", choices = NULL, multiple = TRUE),
                      plotOutput("ProtBoxPlot"),
                      selectInput("level17", 
                                  "Level:", 
                                  choices = c("Protein", "Phosphosite"))
             ),
             
#### Tab18 - Lineplot, single protein ####
             tabPanel("Protein Line",
                      selectizeInput("protein18", "Select Proteins:", choices = NULL, multiple = TRUE),
                      selectizeInput("conditions18", "Select Conditions:", choices = NULL, multiple = TRUE),
                      plotOutput("ProtLinePlot"),
                      selectInput("level18", 
                                  "Level:", 
                                  choices = c("Protein", "Phosphosite"))
             ),
)),
#### SuperTab06 - Phospho-specific #####
tabPanel("Phospho-specific",
         tabsetPanel(
#### Tab19 - Phossite Plot ####
          tabPanel(
            "Phossite Plot",
            sidebarLayout(
              sidebarPanel(
                width = 3,
                numericInput("cutoff19", "Cutoff:", value = 0),
                hr(),
                h4("Annotation Text"),
                selectInput("textPosition19", "Position:", choices = c("Above" = "up", "Below" = "down")),
                textAreaInput("text19", "Enter text:", value = "", width = '100%', height = '100px'),
                actionButton("addText19", "Add"),
                actionButton("deleteText19", "Delete")
              ),
              mainPanel(
                h4("Phossite Plot"),
                customSpinnerWrapper(
                  plotOutput("PhossitePlot", height = "600px"),
                  workmode = workmode
                )
              )
            )
          ),
#### Tab21 - Phossite Coverage Plot ####
          tabPanel(
            "Phossite Coverage Plot",
            sidebarLayout(
              sidebarPanel(
                width = 3,
                actionButton("toggle_id21", "Toggle IDs"),
                hr(),
                h4("Plot Size & Resolution"),
                numericInput("plotWidth21", "Width (cm):", value = 20),
                numericInput("plotHeight21", "Height (cm):", value = 10),
                numericInput("plotDPI21", "DPI:", value = 300),
                downloadButton("downloadPhossitePlot", "Download Plot"),
                hr(),
                h4("Annotation Text"),
                selectInput("textPosition21", "Position:", choices = c("Above" = "up", "Below" = "down")),
                textAreaInput("text21", "Enter text:", value = "", width = '100%', height = '100px'),
                actionButton("addText21", "Add"),
                actionButton("deleteText21", "Delete")
              ),
              mainPanel(
                h4("Phossite Coverage Plot"),
                customSpinnerWrapper(
                  plotOutput("PhossiteCoveragePlot", height = "600px"),
                  workmode = workmode
                )
              )
            )
          ),
#### Tab20.1 - KSEA ####
          tabPanel("KSEA",
                   fluidRow(
                     column(width = 6,
                            selectInput("condition1_20", "Select Condition 1:", choices = c())
                     ),
                     column(width = 6,
                            selectInput("condition2_20", "Select Condition 2:", choices = c())
                     )
                   ),
                   checkboxInput("include_pval_ksea", "Include adjusted p-values in KSEA input", value = FALSE),
                   downloadButton("downKSEA", "Download KSEA input data"),
                   br(), br(),
                   p("1. Select the conditions you want to compare"),
                   p("2. Download the data table"),
                   p("3. Visit the following website:"),
                   tags$a(href = "https://www.phosphosite.org/kinaseLibraryAction", "Website for KSEA", target = "_blank"),
                   br(), br(),
                   p("4. Go to the tab >Fisher Enrichment Analysis<, upload your data and click >Run Enrichment Analysis<"),
                   p("5. Download the table"),
                   fileInput("KSEA_data", "Upload the file here", accept = c(".txt")),
                   plotlyOutput("KinaseVolcPlot"),
                   numericInput("in_pval20", "P-value threshold: ", value = 0.1),
                   downloadButton("downTREE", "Download Kinase Tree input data"),
                   br(), br(),
                   p("6. Visit the following website:"),
                   tags$a(href = "http://phanstiel-lab.med.unc.edu/CORAL/", "Website for Kinase Trees", target = "_blank"),
                   br(), br(),
                   p("7. On the website:"),
                   p("Branch Color > Color Scheme > Select Quantative > Kinases & Value > Paste enrichment_scores.txt values"),
                   p("Node Color > Color Scheme > Select Quantative > Kinases & Value > Paste enrichment_scores.txt values"),
                   p("Node Size > Scaling Scheme > Select Quantative > Kinases & Value > Paste pvals.txt values"),
                   p("Node Size > Missing Kinases > Tick hide"),
                   fileInput("TreeSVG", "Upload and Display Kinase Tree SVG Image", accept = c("image/svg+xml")),
                   uiOutput("treePlot")
          ),

#### Tab20.2 - KSEA (Kinact) ####
              tabPanel("KSEA (Kinact)",
                       fluidRow(
                         column(width = 6,
                                selectInput("condition1_202", "Select Condition 1:", 
                                            choices = c())
                         ),
                         column(width = 6,
                                selectInput("condition2_202", "Select Condition 2:", 
                                            choices = c())
                         )
                       ),
                       plotlyOutput("KinactPlot"),
                       fluidRow(
                         column(width = 6,
                                numericInput("top_n202", "Display Top or Bottom n Kinases: ", value = 0)
                         ),
                         column(width = 6,
                                numericInput("m.cutoff202", "Minimal number of Phosphosites identified in your Dataset for a Kinase: ", value = 5)
                         )
                       ),
                       fluidRow(
                         column(width = 6,
                                selectInput("NetworKIN202", 
                                            "Include NetworKIN prediction?", 
                                            choices = c(FALSE, TRUE))
                         ),
                         column(width = 6,
                                numericInput("NetworKIN.cutoff202", "Minimun NetworKIN score: ", value = 1)
                         )
                       ),
                       textInput("Kinase_202", "Enter a Kinase:"),
                       plotlyOutput("KinactPlotSites"),
                ),

)),
#### SuperTab07 - Tables ####
tabPanel("Tables",
         tabsetPanel(            
#### TabSummary - Summary ####
             tabPanel("Summary", 
                      h4("Data object"),
                      DT::dataTableOutput("transformed_data_prot"),
                      h4("Meta object"),
                      DT::dataTableOutput("meta_object"),
                      downloadButton("download_meta", "Download Metadata"),
                      h4("Phos data object"),
                      DT::dataTableOutput("collapse_data_sum"),
                      downloadButton("download_collapse_data_sum", "Download Collapse Data"),
                      h4("Meta object phos"),
                      DT::dataTableOutput("meta_object2"),
                      downloadButton("download_meta2", "Download Metadata Phos"),
             ),
             
#### Log ####
             #Log
             tabPanel("Log",
                      tableOutput("LogTable"),
                      downloadButton("log", "Download Log"),
                      hr(),
                      fileInput("upload_log", "Upload Log File", accept = ".csv")
             ),
)),
      )
    )
  )
)

#### Server functions ####
server <- function(input, output, session) {
  options(shiny.maxRequestSize=10*1024*1024^2)
#### Tab01 - Data upload ####
  data <- reactive({
    req(input$file)
    
    ext <- tools::file_ext(input$file$name)
    df <- switch(ext,
                 csv = read_csv(input$file$datapath),
                 tsv = read_tsv(input$file$datapath),
                 txt = read_delim(input$file$datapath, delim = "\t"),
                 xlsx = read_excel(input$file$datapath),
                 stop("Invalid file type")
    )
    
    df = rename_cols(df)
    return(df)
  })
  
  data2 <- reactive({
    req(input$file2)
    
    ext <- tools::file_ext(input$file2$name)
    df <- switch(ext,
                 csv = read_csv(input$file2$datapath),
                 tsv = read_tsv(input$file2$datapath),
                 txt = read_delim(input$file2$datapath, delim = "\t"),
                 xlsx = read_excel(input$file2$datapath),
                 stop("Invalid file type")
    )
    
    df = rename_cols(df)
    return(df)
  })
  
  data3 <- reactive({
    req(input$file3)
    
    ext <- tools::file_ext(input$file3$name)
    df <- switch(ext,
                 csv = read_csv(input$file3$datapath),
                 tsv = read_tsv(input$file3$datapath),
                 txt = read_delim(input$file3$datapath, delim = "\t"),
                 xlsx = read_excel(input$file3$datapath),
                 stop("Invalid file type")
    )
    if (input$collapse_option == "not_collapsed"){
      python_script = "Collapse.py"
      command <- sprintf("python %s %s %s", python_script, input$file3$name, input$collapse_cutoff)
      system(command)
      df = read.csv(paste0(getwd(), "/Data/", file_path_sans_ext(input$file3$name), "_collapsed", ".csv"))
    }
    return(df)
  })
  
  output$table <- DT::renderDataTable({
    req(data())
    DT::datatable(data(), options = list(pageLength = 5))
  })
  
  output$table2 <- DT::renderDataTable({
    req(data2())
    DT::datatable(data2(), options = list(pageLength = 5))
  })
  
  output$table3 <- DT::renderDataTable({
    req(data3())
    DT::datatable(data3(), options = list(pageLength = 5))
  })
  
#### Tab02 - Data annotation ####
  #Normal
  output$condition_inputs <- renderUI({
    req(data())
    num_conditions <- input$num_conditions
    lapply(1:num_conditions, function(i) {
      fluidRow(
        column(width = 6,
               textInput(inputId = paste0("condition_", i), label = paste("Condition", i))
        ),
        column(width = 6,
               selectizeInput(
                 inputId = paste0("columns_", i),
                 label = paste("Select columns for Condition", i),
                 choices = colnames(data()),
                 selected = NULL,
                 multiple = TRUE,
                 options = list(
                   dropdownParent = 'body'
                 )
               ),
               tags$style(HTML(sprintf(
                 "#columns_%d.selectize-control { width: 100%% !important; }", i
               )))
        )
      )
    })
  })
  
  observe({
    lapply(1:input$num_conditions, function(i) {
      condition_id <- paste0("condition_", i)
      observeEvent(input[[condition_id]], {
        updateSelectizeInput(session, paste0("columns_", i), label = paste("Select columns for", input[[condition_id]]))
      })
    })
  })
  
  meta <- reactive({
    req(data())
    if (!is.null(input$upload_meta)) {
      file_ext <- tools::file_ext(input$upload_meta$name)
      
      if (file_ext == "csv") {
        meta_df <- read.csv(input$upload_meta$datapath, stringsAsFactors = FALSE)
      } else if (file_ext == "xlsx") {
        meta_df <- read_excel(input$upload_meta$datapath)
      } else if (file_ext == "txt") {
        meta_df <- read_delim(input$upload_meta$datapath, delim = "\t", col_types = cols())
      } else if (file_ext == "tsv") {
        meta_df <- read_data(input$upload_meta$datapath)
      } else {
        stop("Wrong file type!")
      }
      
      if (grepl("Setup", input$upload_meta$name)) {
        meta_df <- transform_meta(data(), meta_df)
      }
      
      return(meta_df)
    }
    
    if (!is.null(input$num_conditions) && input$num_conditions > 0) {
      conditions <- lapply(1:input$num_conditions, function(i) {
        cond_name <- input[[paste0("condition_", i)]]
        cols <- input[[paste0("columns_", i)]]
        if (!is.null(cond_name) && !is.null(cols) && length(cols) > 0) {
          return(data.frame(
            sample = cols,
            condition = cond_name,
            stringsAsFactors = FALSE
          ))
        } else {
          return(NULL)
        }
      })
      conditions <- do.call(rbind, conditions)
      if (!is.null(conditions)) return(conditions)
    }
    
    sample_columns <- colnames(data())[grepl("\\[", colnames(data()))]
    if (length(sample_columns) > 0) {
      return(data.frame(
        sample = sample_columns,
        condition = "sample",
        stringsAsFactors = FALSE
      ))
    }
    
    return(NULL)
  })
  
  
  was_transformed <- reactiveVal(FALSE)
  transformed_data <- reactiveVal(NULL)
  not_transformed_data <- reactiveVal(NULL)
  
  observeEvent(input$log2_yes, {
    tryCatch({
      req(data(), meta())
      transformed_data(data())
      not_transformed_data(inverseof_log2_transform_data(data(), meta()))
      was_transformed(FALSE)
    }, error = function(e) {
      showNotification("Meta data is missing or invalid. Please define it before applying log2 settings.", type = "error")
    })
  })
  
  observeEvent(input$log2_no, {
    tryCatch({
      req(data(), meta())
      transformed_data(log2_transform_data(data(), meta()))
      not_transformed_data(data())
      was_transformed(TRUE)
    }, error = function(e) {
      showNotification("Meta data is missing or invalid. Please define it before applying log2 settings.", type = "error")
    })
  })
  
  
  filtered_data <- reactiveVal(NULL)
  log2_filtered_data <- reactiveVal(NULL)
  
  observeEvent(input$apply_filter, {
    req(data(), transformed_data(), meta(), input$filter_num)
    filtered_data(filter_data(data(), meta(), num = input$filter_num, filterops = input$filterop1))
    log2_filtered_data(filter_data(transformed_data(), meta(), num = input$filter_num, filterops = input$filterop1))
  })
  
  
  output$displayed_data <- DT::renderDataTable({
    if (is.null(log2_filtered_data())) {
      req(transformed_data())
      DT::datatable(transformed_data(), options = list(pageLength = 5))
    } else {
      req(log2_filtered_data())
      DT::datatable(log2_filtered_data(), options = list(pageLength = 5))
    }
  })
  
  pre_volcano_data <- reactive({
      if (!is.null(log2_filtered_data())) {
        return(log2_filtered_data())
      } else {
        req(transformed_data())
        return(transformed_data())
      }
    })
  
  volcano_data <- reactive({
    if (!is.null(imputed_data())) {
      req(imputed_data())
      return(imputed_data())
    } else {
      if (!is.null(log2_filtered_data())) {
        return(log2_filtered_data())
      } else {
        req(transformed_data())
        return(transformed_data())
      }
    }
  })
  
  #Phospho
  output$condition_inputs2 <- renderUI({
    req(data3())
    num_conditions2 <- input$num_conditions2
    lapply(1:num_conditions2, function(i) {
      fluidRow(
        column(width = 6,
               textInput(inputId = paste0("condition2_", i), label = paste("Condition", i))
        ),
        column(width = 6,
               selectizeInput(
                 inputId = paste0("columns2_", i),
                 label = paste("Select columns for Condition", i),
                 choices = colnames(data3()),
                 selected = NULL,
                 multiple = TRUE,
                 options = list(
                   dropdownParent = 'body'
                 )
               ),
               tags$style(HTML(sprintf(
                 "#columns2_%d.selectize-control { width: 100%% !important; }", i
               )))
        )
      )
    })
  })
  
  observe({
    lapply(1:input$num_conditions2, function(i) {
      condition_id2 <- paste0("condition2_", i)
      observeEvent(input[[condition_id2]], {
        updateSelectizeInput(session, paste0("columns2_", i), label = paste("Select columns for", input[[condition_id2]]))
      })
    })
  })
  
  meta2 <- reactive({
    req(data3())
    if (!is.null(input$upload_meta2)) {
      file_ext <- tools::file_ext(input$upload_meta2$name)
      
      if (file_ext == "csv") {
        meta_df2 <- read.csv(input$upload_meta2$datapath, stringsAsFactors = FALSE)
      } else if (file_ext == "xlsx") {
        meta_df2 <- read_excel(input$upload_meta2$datapath)
      } else if (file_ext == "txt") {
        meta_df2 <- read_delim(input$upload_meta2$datapath, delim = "\t", col_types = cols())
      } else if (file_ext == "tsv") {
        meta_df2 <- read_data(input$upload_meta2$datapath)
      } else {
        stop("Wrong file type!")
      }
      
      if (grepl("Setup", input$upload_meta2$name)) {
        meta_df2 <- transform_meta(data3(), meta_df2)
      }
      
      return(meta_df2)
    }
    
    if (!is.null(input$num_conditions2) && input$num_conditions2 > 0) {
      conditions2 <- lapply(1:input$num_conditions2, function(i) {
        cond_name <- input[[paste0("condition2_", i)]]
        cols <- input[[paste0("columns2_", i)]]
        if (!is.null(cond_name) && !is.null(cols) && length(cols) > 0) {
          return(data.frame(
            sample = cols,
            condition = cond_name,
            stringsAsFactors = FALSE
          ))
        } else {
          return(NULL)
        }
      })
      conditions2 <- do.call(rbind, conditions2)
      if (!is.null(conditions2)) return(conditions2)
    }
    
    sample_columns2 <- colnames(data3())[grepl("\\[", colnames(data3()))]
    if (length(sample_columns2) > 0) {
      return(data.frame(
        sample = sample_columns2,
        condition = "sample",
        stringsAsFactors = FALSE
      ))
    }
    
    return(NULL)
  })
  
  
  was_transformed3 <- reactiveVal(FALSE)
  transformed_data3 <- reactiveVal(NULL)
  not_transformed_data3 <- reactiveVal(NULL)
  
  observeEvent(input$log2_yes2, {
    tryCatch({
      req(data3(), meta2())
      transformed_data3(data3())
      not_transformed_data3(inverseof_log2_transform_data(data3(), meta2()))
      was_transformed3(FALSE)
    }, error = function(e) {
      showNotification("Meta data is missing or invalid. Please define it before applying log2 settings.", type = "error")
    })
  })
  
  observeEvent(input$log2_no2, {
    tryCatch({
      req(data3(), meta2())
      transformed_data3(log2_transform_data(data3(), meta2()))
      not_transformed_data3(data3())
      was_transformed3(TRUE)
    }, error = function(e) {
      showNotification("Meta data is missing or invalid. Please define it before applying log2 settings.", type = "error")
    })
  })
  
  filtered_data3 <- reactiveVal(NULL)
  log2_filtered_data3 <- reactiveVal(NULL)

  observeEvent(input$apply_filter2, {
    req(data3(), transformed_data3(), meta2(), input$filter_num2)
    filtered_data3(filter_data(data3(), meta2(), num = input$filter_num2, filterops = input$filterop2))
    log2_filtered_data3(filter_data(transformed_data3(), meta2(), num = input$filter_num2, filterops = input$filterop2))
  })  
  
  output$displayed_data2 <- DT::renderDataTable({
    if (is.null(log2_filtered_data3())) {
      req(transformed_data3())
      DT::datatable(transformed_data3(), options = list(pageLength = 5))
    } else {
      req(log2_filtered_data3())
      DT::datatable(log2_filtered_data3(), options = list(pageLength = 5))
    }
  })
  
  pre_volcano_data3 <- reactive({
    if (!is.null(log2_filtered_data3())) {
      return(log2_filtered_data3())
    } else {
      req(transformed_data3())
      return(transformed_data3())
    }
  })
  
  volcano_data3 <- reactive({
    if (!is.null(imputed_data3())) {
      req(imputed_data3())
      return(imputed_data3())
    } else {
      if (!is.null(log2_filtered_data3())) {
        return(log2_filtered_data3())
      } else {
        req(transformed_data3())
        return(transformed_data3())
      }
    }
  })
  
  #Color
  observeEvent(input$color_palette, {
    init_colors(color = input$color_palette)
  })
  
  observeEvent(input$reloadButton, {
    #CoveragePlot
    output$coveragePlot <- renderPlot({
      if (input$level3 == "Protein") {
        req(data(), meta())
        coverage_plot(data(), meta(), id3())
      } else if (input$level3 == "Peptide") {  
        req(data2(), meta())
        coverage_plot_pep(data2(), meta(), id3())
      } else if (input$level3 == "Phosphosite") {  
        req(data3(), meta2())
        coverage_plot(data3(), meta2(), id3())
      }
    })
    #HistoInt
    output$HistIntPlot <- renderPlot({
      if (input$level5 == "Protein") {
        req(transformed_data(), meta())
        histo_int(transformed_data(), meta())
      } else if (input$level5 == "Phosphosite"){
        req(data3(), meta2())
        histo_int(data3(), meta2())
      }
    })
    #BoxplotInt
    output$BoxIntPlot <- renderPlot({
      if (input$level6 == "Protein") {
        req(transformed_data(), meta())
        if (meanBX()){
          boxplot_int(transformed_data(), meta(), outliers = outliersBX())
        }else{
          boxplot_int_single(transformed_data(), meta(), outliers = outliersBX(), id=id6()) 
        }
      } else if (input$level6 == "Phosphosite"){
        req(data3(), meta2())
        if (meanBX()){
          boxplot_int(data3(), meta2(), outliers = outliersBX())
        }else{
          boxplot_int_single(data3(), meta2(), outliers = outliersBX(), id=id6()) 
        }
      }
    })
    #Cov Plot
    output$CovPlot <- renderPlot({
      if (input$level7 == "Protein") {
        req(data(), meta())
        cov_plot(not_transformed_data(), meta(), outliers = outliersCOV())
      } else if (input$level7 == "Phosphosite"){
        req(data3(), meta2())
        cov_plot(not_transformed_data3(), meta2(), outliers = outliersCOV())
      }
    })
    #PCA Plot
    output$PCAPlot <- renderPlot({
      if (input$level8 == "Protein") {
        req(data(), meta())
        dim_func(data(), meta(), method=input$style8)
      } else if (input$level8 == "Phosphosite"){
        req(data3(), meta2())
        dim_func(data3(), meta2(), method=input$style8)
      }
    })
    #Abundance Plot
    output$abundancePlot <- renderPlotly({
      if (input$level9 == "Protein"){
        req(input$condition9)
        if (input$condition9 == "All Conditions") {
          abundance_plot(not_transformed_data(), meta(), workflow = input$level9)
        } else {
          interactive_abundance_plot(not_transformed_data(), meta(), input$condition9, workflow = input$level9, search = input$protein9)
        }} else if (input$level9 == "Phosphosite"){
          req(input$condition9)
          if (input$condition9 == "All Conditions") {
            abundance_plot(not_transformed_data3(), meta2(), workflow = "Phosphosite")
          } else {
            interactive_abundance_plot(not_transformed_data3(), meta2(), input$condition9, workflow = "Phosphosite", search = input$protein9)
          }}
    })
  })
  
#### Tab3.5 - Data Distribution ####
  output$FirstDigitPlot <- renderPlot({
    if (input$level3.5 == "Protein") {
      req(not_transformed_data(), meta())
      first_digit_distribution(not_transformed_data(), meta())
    } else if (input$level3.5 == "Phosphosite") {
      req(data3(), meta2())
      first_digit_distribution(not_transformed_data3(), meta2())
    }
  })
  
  output$DataStructurePlot <- renderPlot({
    if (input$level3.5 == "Protein") {
      req(not_transformed_data(), meta())
      data_pattern_structure(not_transformed_data(), meta())
    } else if (input$level3.5 == "Phosphosite") {
      req(data3(), meta2())
      data_pattern_structure(not_transformed_data3(), meta2())
    }
  })
  
#### Tab4.5 - Data Distribution ####
  output$qqnorm <- renderPlot({
    if (input$level4.5 == "Protein") {
      req(transformed_data(), meta())
      qqnorm_plot(transformed_data(), meta())
    } else if (input$level4.5 == "Phosphosite") {
      req(data3(), meta2())
      qqnorm_plot(transformed_data3(), meta2())
    }
  })
  
#### Tab03 - Coverage Plot ####
  id3 <- reactiveVal(TRUE)
  header3 <- reactiveVal(TRUE)
  legend3 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_id3, {
    id3(!id3())
  })
  
  observeEvent(input$toggle_header3, {
    header3(!header3())
  })
  
  observeEvent(input$toggle_legend3, {
    legend3(!legend3())
  })
  
  output$coveragePlot <- renderPlot({
    if (input$level3 == "Protein") {
      req(data(), meta())
      coverage_plot(data(), meta(), id = id3(), header = header3(), legend = legend3())
    } else if (input$level3 == "Peptide") {
      req(data2(), meta())
      coverage_plot_pep(data2(), meta(), id = id3(), header = header3(), legend = legend3())
    } else if (input$level3 == "Phosphosite") {
      req(data3(), meta2())
      coverage_plot(data3(), meta2(), id = id3(), header = header3(), legend = legend3())
    }
  })
  
  output$downloadCoveragePlot <- downloadHandler(
    filename = function() {
      ext <- input$plotFormat3
      paste("coverage_plot", Sys.Date(), ".", ext, sep = "")
    },
    content = function(file) {
      width <- input$plotWidth3
      height <- input$plotHeight3
      dpi <- input$plotDPI3
      ext <- input$plotFormat3
      
      p <- NULL
      if (input$level3 == "Protein") {
        req(data(), meta())
        p <- coverage_plot(data(), meta(), id = id3(), header = header3(), legend = legend3())
      } else if (input$level3 == "Peptide") {
        req(data2(), meta())
        p <- coverage_plot_pep(data2(), meta(), id = id3(), header = header3(), legend = legend3())
      } else if (input$level3 == "Phosphosite") {
        req(data3(), meta2())
        p <- coverage_plot(data3(), meta2(), id = id3(), header = header3(), legend = legend3())
      }
      
      if (ext %in% c("png", "jpg")) {
        ggsave(file, plot = p, device = ext, width = width, height = height, units = "cm", dpi = dpi)
      } else if (ext == "pdf") {
        pdf(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      } else if (ext == "svg") {
        svglite::svglite(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      }
    }
  )
  
  observeEvent(input$addText3, {
    df <- const_df()
    suffix <- if (input$textPosition3 == "up") "text3up" else "text3down"
    varName <- switch(input$level3,
                      "Protein" = paste0(suffix, "Protein"),
                      "Peptide" = paste0(suffix, "Peptide"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"))
    df <- df[df$Var != varName, ]
    new_row <- data.frame(Var = varName, Select = input$text3, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text3", value = "")
  })
  
  observeEvent(input$deleteText3, {
    df <- const_df()
    suffix <- if (input$textPosition3 == "up") "text3up" else "text3down"
    varName <- switch(input$level3,
                      "Protein" = paste0(suffix, "Protein"),
                      "Peptide" = paste0(suffix, "Peptide"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"))
    df <- df[df$Var != varName, ]
    const_df(df)
  })
  
  observe({
    df <- const_df()
    suffix <- if (input$textPosition3 == "up") "text3up" else "text3down"
    varName <- switch(input$level3,
                      "Protein" = paste0(suffix, "Protein"),
                      "Peptide" = paste0(suffix, "Peptide"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"))
    
    value <- df$Select[df$Var == varName]
    updateTextAreaInput(session, "text3", value = if (length(value) > 0) value else "")
  })
  
#### Tab04 - Missing Value Plot ####
  header4 <- reactiveVal(TRUE)
  text4 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_header4, {
    header4(!header4())
  })
  
  observeEvent(input$toggle_text4, {
    text4(!text4())
  })
  
  output$MissValPlot <- renderPlot({
    bin4 <- input$missValBin4
    
    p <- NULL
    if (input$level4 == "Protein") {
      req(data(), meta())
      p <- missing_value_plot(data=data(), meta=meta(), bin=bin4, text=text4(), text_size=input$text_size4)
    } else if (input$level4 == "Precursor") {
      req(data2(), meta())
      p <- missing_value_plot_prec(data2=data2(), meta=meta(), bin=bin4, text=text4(), text_size=input$text_size4)
    } else if (input$level4 == "Peptide") {
      req(data2(), meta())
      p <- missing_value_plot_pep(data2=data2(), meta=meta(), bin=bin4, text=text4(), text_size=input$text_size4)
    } else if (input$level4 == "Phosphosite") {
      req(data3(), meta2())
      p <- missing_value_plot(data=data3(), meta=meta2(), bin=bin4, text=text4(), text_size=input$text_size4)
    }
    
    if (!header4()) {
      p <- p + ggtitle(NULL)
    }
    
    p
  })
  
  output$downloadMissValPlot <- downloadHandler(
    filename = function() {
      ext <- input$missValPlotFormat
      paste0("missval_plot_", Sys.Date(), ".", ext)
    },
    content = function(file) {
      width <- input$plotWidth4
      height <- input$plotHeight4
      dpi <- input$plotDPI4
      bin4 <- input$missValBin4
      ext <- input$missValPlotFormat
      
      plot_to_save <- switch(input$level4,
                             "Protein" = missing_value_plot(data(), meta(), bin4),
                             "Precursor" = missing_value_plot_prec(data2(), meta(), bin4),
                             "Peptide" = missing_value_plot_pep(data2(), meta(), bin4),
                             "Phosphosite" = missing_value_plot(data3(), meta2(), bin4)
      )
      
      if (!header4()) {
        plot_to_save <- plot_to_save + ggtitle(NULL)
      }
      
      ggsave(
        filename = file,
        plot = plot_to_save,
        device = ext,
        width = width,
        height = height,
        units = "cm",
        dpi = dpi
      )
    }
  )
  
  observeEvent(input$addText4, {
    df <- const_df()
    suffix <- if (input$textPosition4 == "up") "text4up" else "text4down"
    varName <- switch(input$level4,
                      "Protein" = paste0(suffix, "Protein"),
                      "Peptide" = paste0(suffix, "Peptide"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"),
                      "Precursor" = paste0(suffix, "Precursor"))
    df <- df[df$Var != varName, ]
    new_row <- data.frame(Var = varName, Select = input$text4, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text4", value = "")
  })
  
  observeEvent(input$deleteText4, {
    df <- const_df()
    suffix <- if (input$textPosition4 == "up") "text4up" else "text4down"
    varName <- switch(input$level4,
                      "Protein" = paste0(suffix, "Protein"),
                      "Peptide" = paste0(suffix, "Peptide"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"),
                      "Precursor" = paste0(suffix, "Precursor"))
    df <- df[df$Var != varName, ]
    const_df(df)
  })
  
  observe({
    df <- const_df()
    suffix <- if (input$textPosition4 == "up") "text4up" else "text4down"
    varName <- switch(input$level4,
                      "Protein" = paste0(suffix, "Protein"),
                      "Peptide" = paste0(suffix, "Peptide"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"),
                      "Precursor" = paste0(suffix, "Precursor"))
    value <- df$Select[df$Var == varName]
    updateTextAreaInput(session, "text4", value = if (length(value) > 0) value else "")
  })
  
#### Tab05 - Histogram ####
  header5 <- reactiveVal(TRUE)
  legend5 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_header5, {
    header5(!header5())
  })
  
  observeEvent(input$toggle_legend5, {
    legend5(!legend5())
  })
  
  output$HistIntPlot <- renderPlot({
    if (input$level5 == "Protein") {
      req(transformed_data(), meta())
      histo_int(transformed_data(), meta(), header = header5(), legend = legend5())
    } else if (input$level5 == "Phosphosite") {
      req(data3(), meta2())
      histo_int(data3(), meta2(), header = header5(), legend = legend5())
    }
  })
  
  output$downloadHistIntPlot <- downloadHandler(
    filename = function() {
      ext <- input$plotFormat5
      paste("histint_plot", Sys.Date(), ".", ext, sep = "")
    },
    content = function(file) {
      width <- input$plotWidth5
      height <- input$plotHeight5
      dpi <- input$plotDPI5
      ext <- input$plotFormat5
      
      p <- NULL
      if (input$level5 == "Protein") {
        req(transformed_data(), meta())
        p <- histo_int(transformed_data(), meta(), header = header5(), legend = legend5())
      } else if (input$level5 == "Phosphosite") {
        req(data3(), meta2())
        p <- histo_int(data3(), meta2(), header = header5(), legend = legend5())
      }
      
      if (ext %in% c("png", "jpg")) {
        ggsave(file, plot = p, device = ext, width = width, height = height, units = "cm", dpi = dpi)
      } else if (ext == "pdf") {
        pdf(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      } else if (ext == "svg") {
        svglite::svglite(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      }
    }
  )
  
  observeEvent(input$addText5, {
    df <- const_df()
    suffix <- if (input$textPosition5 == "up") "text5up" else "text5down"
    varName <- switch(input$level5,
                      "Protein" = paste0(suffix, "Protein"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"))
    df <- df[df$Var != varName, ]
    new_row <- data.frame(Var = varName, Select = input$text5, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text5", value = "")
  })
  
  observeEvent(input$deleteText5, {
    df <- const_df()
    suffix <- if (input$textPosition5 == "up") "text5up" else "text5down"
    varName <- switch(input$level5,
                      "Protein" = paste0(suffix, "Protein"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"))
    df <- df[df$Var != varName, ]
    const_df(df)
  })
  
  observe({
    df <- const_df()
    suffix <- if (input$textPosition5 == "up") "text5up" else "text5down"
    varName <- switch(input$level5,
                      "Protein" = paste0(suffix, "Protein"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"))
    value <- df$Select[df$Var == varName]
    updateTextAreaInput(session, "text5", value = if (length(value) > 0) value else "")
  })
  
#### Tab06 - Boxplot ####
  id6 <- reactiveVal(TRUE)
  outliersBX <- reactiveVal(FALSE)
  meanBX <- reactiveVal(FALSE)
  header6 <- reactiveVal(TRUE)
  legend6 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_id6, {
    id6(!id6())
  })
  
  observeEvent(input$toggle_outliersBX, {
    outliersBX(!outliersBX())
  })
  
  observeEvent(input$toggle_meanBX, {
    meanBX(!meanBX())
  })
  
  observeEvent(input$toggle_header6, {
    header6(!header6())
  })
  
  observeEvent(input$toggle_legend6, {
    legend6(!legend6())
  })
  
  output$BoxIntPlot <- renderPlot({
    if (input$level6 == "Protein") {
      req(transformed_data(), meta())
      if (meanBX()) {
        boxplot_int(transformed_data(), meta(), outliers = outliersBX(), header = header6(), legend = legend6())
      } else {
        boxplot_int_single(transformed_data(), meta(), outliers = outliersBX(), id = id6(), header = header6(), legend = legend6())
      }
    } else if (input$level6 == "Phosphosite") {
      req(data3(), meta2())
      if (meanBX()) {
        boxplot_int(data3(), meta2(), outliers = outliersBX(), header = header6(), legend = legend6())
      } else {
        boxplot_int_single(data3(), meta2(), outliers = outliersBX(), id = id6(), header = header6(), legend = legend6())
      }
    }
  })
  
  output$downloadBoxIntPlot <- downloadHandler(
    filename = function() {
      ext <- input$plotFormat6
      paste("boxint_plot", Sys.Date(), ".", ext, sep = "")
    },
    content = function(file) {
      width <- input$plotWidth6
      height <- input$plotHeight6
      dpi <- input$plotDPI6
      ext <- input$plotFormat6
      
      p <- if (input$level6 == "Protein") {
        req(transformed_data(), meta())
        if (meanBX()) {
          boxplot_int(transformed_data(), meta(), outliers = outliersBX(), header = header6(), legend = legend6())
        } else {
          boxplot_int_single(transformed_data(), meta(), outliers = outliersBX(), id = id6(), header = header6(), legend = legend6())
        }
      } else {
        req(data3(), meta2())
        if (meanBX()) {
          boxplot_int(data3(), meta2(), outliers = outliersBX(), header = header6(), legend = legend6())
        } else {
          boxplot_int_single(data3(), meta2(), outliers = outliersBX(), id = id6(), header = header6(), legend = legend6())
        }
      }
      
      if (ext %in% c("png", "jpg")) {
        ggsave(file, plot = p, device = ext, width = width, height = height, units = "cm", dpi = dpi)
      } else if (ext == "pdf") {
        pdf(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      } else if (ext == "svg") {
        svglite::svglite(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      }
    }
  )
  
  observeEvent(input$addText6, {
    df <- const_df()
    suffix <- if (input$textPosition6 == "up") "text6up" else "text6down"
    varName <- switch(input$level6,
                      "Protein" = paste0(suffix, "Protein"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"))
    df <- df[df$Var != varName, ]
    new_row <- data.frame(Var = varName, Select = input$text6, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text6", value = "")
  })
  
  observeEvent(input$deleteText6, {
    df <- const_df()
    suffix <- if (input$textPosition6 == "up") "text6up" else "text6down"
    varName <- switch(input$level6,
                      "Protein" = paste0(suffix, "Protein"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"))
    df <- df[df$Var != varName, ]
    const_df(df)
  })
  
  observe({
    df <- const_df()
    suffix <- if (input$textPosition6 == "up") "text6up" else "text6down"
    varName <- switch(input$level6,
                      "Protein" = paste0(suffix, "Protein"),
                      "Phosphosite" = paste0(suffix, "Phosphosite"))
    value <- df$Select[df$Var == varName]
    updateTextAreaInput(session, "text6", value = if (length(value) > 0) value else "")
  })
  
#### Tab07 - COV Plot ####
  outliersCOV <- reactiveVal(FALSE)
  header7 <- reactiveVal(TRUE)
  legend7 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_outliersCOV, {
    outliersCOV(!outliersCOV())
  })
  
  observeEvent(input$toggle_header7, {
    header7(!header7())
  })
  
  observeEvent(input$toggle_legend7, {
    legend7(!legend7())
  })
  
  output$CovPlot <- renderPlot({
    if (input$level7 == "Protein") {
      req(data(), meta())
      cov_plot(not_transformed_data(), meta(), outliers = outliersCOV(), header = header7(), legend = legend7())
    } else if (input$level7 == "Phosphosite") {
      req(data3(), meta2())
      cov_plot(not_transformed_data3(), meta2(), outliers = outliersCOV(), header = header7(), legend = legend7())
    }
  })
  
  output$downloadCovPlot <- downloadHandler(
    filename = function() {
      paste0("cov_plot_", Sys.Date(), ".", input$plotFormat7)
    },
    content = function(file) {
      width <- input$plotWidth7
      height <- input$plotHeight7
      dpi <- input$plotDPI7
      ext <- input$plotFormat7
      
      p <- if (input$level7 == "Protein") {
        req(not_transformed_data(), meta())
        cov_plot(not_transformed_data(), meta(), outliers = outliersCOV(), header = header7(), legend = legend7())
      } else {
        req(not_transformed_data3(), meta2())
        cov_plot(not_transformed_data3(), meta2(), outliers = outliersCOV(), header = header7(), legend = legend7())
      }
      
      if (ext %in% c("png", "jpg")) {
        ggsave(file, plot = p, device = ext, width = width, height = height, units = "cm", dpi = dpi)
      } else if (ext == "pdf") {
        pdf(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      } else if (ext == "svg") {
        svglite::svglite(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      }
    }
  )
  
  observeEvent(input$addText7, {
    df <- const_df()
    suffix <- if (input$textPosition7 == "up") "text7up" else "text7down"
    varName <- if (input$level7 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    df <- df[df$Var != varName, ]
    new_row <- data.frame(Var = varName, Select = input$text7, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text7", value = "")
  })
  
  observeEvent(input$deleteText7, {
    df <- const_df()
    suffix <- if (input$textPosition7 == "up") "text7up" else "text7down"
    varName <- if (input$level7 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    df <- df[df$Var != varName, ]
    const_df(df)
  })
  
  observe({
    df <- const_df()
    suffix <- if (input$textPosition7 == "up") "text7up" else "text7down"
    varName <- if (input$level7 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    value <- df$Select[df$Var == varName]
    updateTextAreaInput(session, "text7", value = if (length(value) > 0) value else "")
  })
  
#### Tab08 - PCA ####
  header8 <- reactiveVal(TRUE)
  legend8 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_header8, {
    header8(!header8())
  })
  
  observeEvent(input$toggle_legend8, {
    legend8(!legend8())
  })
  
  output$PCAPlot <- renderPlot({
    req(input$type8 == "Non-Interactive")
    
    dot_size <- input$dotSize8
    
    if (input$level8 == "Protein") {
      req(data(), meta())
      dim_func(data(), meta(), method = input$style8, header = header8(), legend = legend8(), dot_size = dot_size)
    } else if (input$level8 == "Phosphosite") {
      req(data3(), meta2())
      dim_func(data3(), meta2(), method = input$style8, header = header8(), legend = legend8(), dot_size = dot_size)
    }
  })
  
  output$PCAPlotInteractive <- renderPlotly({
    req(input$type8 == "Interactive")
    
    if (input$level8 == "Protein") {
      req(data(), meta())
      pca_plot_interactive(data(), meta(), header = header8(), legend = legend8())
    } else if (input$level8 == "Phosphosite") {
      req(data3(), meta2())
      pca_plot_interactive(data3(), meta2(), header = header8(), legend = legend8())
    }
  })
  
  output$downloadPCAPlot <- downloadHandler(
    filename = function() {
      paste0("pca_plot_", Sys.Date(), ".", input$plotFormat8)
    },
    content = function(file) {
      width <- input$plotWidth8
      height <- input$plotHeight8
      dpi <- input$plotDPI8
      ext <- input$plotFormat8
      dot_size <- input$dotSize8
      
      p <- NULL
      if (input$level8 == "Protein") {
        req(data(), meta())
        p <- dim_func(data(), meta(), method = input$style8, header = header8(), legend = legend8(), dot_size = dot_size)
      } else if (input$level8 == "Phosphosite") {
        req(data3(), meta2())
        p <- dim_func(data3(), meta2(), method = input$style8, header = header8(), legend = legend8(), dot_size = dot_size)
      }
      
      if (ext %in% c("png", "jpg")) {
        ggsave(file, plot = p, device = ext, width = width, height = height, units = "cm", dpi = dpi)
      } else if (ext == "pdf") {
        pdf(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      } else if (ext == "svg") {
        svglite::svglite(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      }
    }
  )
  
  observeEvent(input$addText8, {
    df <- const_df()
    suffix <- if (input$textPosition8 == "up") "text8up" else "text8down"
    varName <- if (input$level8 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    df <- df[df$Var != varName, ]
    new_row <- data.frame(Var = varName, Select = input$text8, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text8", value = "")
  })
  
  observeEvent(input$deleteText8, {
    df <- const_df()
    suffix <- if (input$textPosition8 == "up") "text8up" else "text8down"
    varName <- if (input$level8 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    df <- df[df$Var != varName, ]
    const_df(df)
  })
  
  observe({
    df <- const_df()
    suffix <- if (input$textPosition8 == "up") "text8up" else "text8down"
    varName <- if (input$level8 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    value <- df$Select[df$Var == varName]
    updateTextAreaInput(session, "text8", value = if (length(value) > 0) value else "")
  })
  
#### Tab09 - Abundance Plot ####
  observe({
    if (input$level9 == "Protein"){
      req(not_transformed_data())
      updateSelectizeInput(session, "protein9", choices = unique(not_transformed_data()$ProteinNames))
    } else if (input$level9 == "Phosphosite"){
      req(not_transformed_data3())
      updateSelectizeInput(session, "protein9", choices = unique(not_transformed_data3()$PTM_Collapse_key))
    }
  })
  
  observe({
    if (input$level9 == "Protein"){
      req(not_transformed_data())
      updateSelectInput(session, "condition9", 
                        choices = c("All Conditions", unique(meta()$condition)))}
    else if (input$level9 == "Phosphosite"){
      req(data3())
      updateSelectInput(session, "condition9", 
                        choices = c("All Conditions", unique(meta2()$condition)))}
  })
  
  output$abundancePlot <- renderPlotly({
    if (input$level9 == "Protein"){
      req(input$condition9)
      if (input$condition9 == "All Conditions") {
        abundance_plot(not_transformed_data(), meta(), workflow = input$level9)
      } else {
        interactive_abundance_plot(not_transformed_data(), meta(), input$condition9, workflow = input$level9, search = input$protein9)
      }} else if (input$level9 == "Phosphosite"){
        req(input$condition9)
        if (input$condition9 == "All Conditions") {
          abundance_plot(not_transformed_data3(), meta2(), workflow = "Phosphosite")
        } else {
          interactive_abundance_plot(not_transformed_data3(), meta2(), input$condition9, workflow = "Phosphosite", search = input$protein9)
      }}
  })
  
  output$downloadAbPlot <- downloadHandler(
    filename = function() {
      paste("abundance_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth9
      height <- input$plotHeight9
      dpi <- input$plotDPI9
      
      if (input$condition9 == "All Conditions") {
        ggsave(file, plot = abundance_plot(data(), meta()), device = "png", 
               width = width, height = height, units = "cm", dpi = dpi)
      } else {
        print("Please work on this!")
      }
    }
  )
  
  observeEvent(input$addText9, {
    df <- const_df()
    suffix <- if (input$textPosition9 == "up") "text9up" else "text9down"
    varName <- if (input$level9 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    df <- df[df$Var != varName, ]
    new_row <- data.frame(Var = varName, Select = input$text9, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text9", value = "")
  })
  
  observeEvent(input$deleteText9, {
    df <- const_df()
    suffix <- if (input$textPosition9 == "up") "text9up" else "text9down"
    varName <- if (input$level9 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    df <- df[df$Var != varName, ]
    const_df(df)
  })
  
  observe({
    df <- const_df()
    suffix <- if (input$textPosition9 == "up") "text9up" else "text9down"
    varName <- if (input$level9 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    value <- df$Select[df$Var == varName]
    updateTextAreaInput(session, "text9", value = if (length(value) > 0) value else "")
  })
  
#### Tab10 - Volcano Plot ####
  observe({
    if (input$level10 == "Protein"){
      req(transformed_data())
      updateSelectInput(session, "condition1_10", 
                        choices = c(unique(meta()$condition)))
    } else if (input$level10=="Phosphosite"){
      req(transformed_data3())
      updateSelectInput(session, "condition1_10", 
                        choices = c(unique(meta2()$condition)))
    }
  })
  
  observe({
    if (input$level10 == "Protein"){
      req(transformed_data())
      updateSelectInput(session, "condition2_10", 
                        choices = c(unique(meta()$condition)))
    } else if (input$level10=="Phosphosite"){
      req(transformed_data3())
      updateSelectInput(session, "condition2_10", 
                        choices = c(unique(meta2()$condition)))
    }
  })
  
  output$VolcPlot <- renderPlotly({
    if (input$level10 == "Protein"){
      req(volcano_data(), meta(), input$condition1_10, input$condition2_10)
      volcano_plot(volcano_data(), meta(), input$condition1_10, input$condition2_10, in_pval = input$in_pval10, in_log2fc = input$in_log2fc10, paired=input$paired10, uncorrected = input$uncorrected10)
    } else if (input$level10 == "Phosphosite"){
      output$VolcPlot <- renderPlotly({
        req(volcano_data3(), meta2(), input$condition1_10, input$condition2_10)
        volcano_plot(volcano_data3(), meta2(), input$condition1_10, input$condition2_10, in_pval = input$in_pval10, in_log2fc = input$in_log2fc10, workflow="Phosphosite", paired=input$paired10, uncorrected = input$uncorrected10)})
    }
  })
  
  data10 <- reactive({
    if (input$level10 == "Protein") {
      req(volcano_data(), meta(), input$condition1_10, input$condition2_10)
      data <- volcano_data_f(volcano_data(), meta(), input$condition1_10, input$condition2_10, in_pval = input$in_pval10, in_log2fc = input$in_log2fc10, paired=input$paired10, uncorrected = input$uncorrected10)
    } else if (input$level10 == "Phosphosite") {
      req(volcano_data3(), meta2(), input$condition1_10, input$condition2_10)
      data <- volcano_data_f(volcano_data3(), meta2(), input$condition1_10, input$condition2_10, in_pval = input$in_pval10, in_log2fc = input$in_log2fc10, workflow="Phosphosite", paired=input$paired10, uncorrected = input$uncorrected10)
    }
    rownames(data) <- NULL
    return(data)
  })
  
  observeEvent(input$addVolc, {
    if (input$level10 == "Protein") {
      req(volcano_data(), meta(), input$condition1_10, input$condition2_10)
      df <- const_df()
      volc_in <- c(
        input$condition1_10,
        input$condition2_10,
        input$in_pval10,
        input$in_log2fc10,
        as.character(input$uncorrected10)
      )
      volc_in <- paste(volc_in, collapse = "ඞ")
      new_row <- data.frame(Var = "VolcPlot", Select = volc_in, stringsAsFactors = FALSE)
      const_df(rbind(df, new_row))
    } else if (input$level10 == "Phosphosite") {
      req(volcano_data3(), meta2(), input$condition1_10, input$condition2_10)
      df <- const_df()
      phosvolc_in <- c(
        input$condition1_10,
        input$condition2_10,
        input$in_pval10,
        input$in_log2fc10,
        as.character(input$uncorrected10)
      )
      phosvolc_in <- paste(phosvolc_in, collapse = "ඞ")
      new_row <- data.frame(Var = "PhosVolcPlot", Select = phosvolc_in, stringsAsFactors = FALSE)
      const_df(rbind(df, new_row))
    }
  })
  
  observeEvent(input$removeVolc, {
    req(input$condition1_10, input$condition2_10)
    df <- const_df()
    current_key <- paste(
      input$condition1_10,
      input$condition2_10,
      input$in_pval10,
      input$in_log2fc10,
      as.character(input$uncorrected10),
      sep = "ඞ"
    )
    if (input$level10 == "Protein") {
      df <- df[!(df$Var == "VolcPlot" & df$Select == current_key), ]
    } else if (input$level10 == "Phosphosite") {
      df <- df[!(df$Var == "PhosVolcPlot" & df$Select == current_key), ]
    }
    const_df(df)
  })
  
  observeEvent(input$addText10, {
    df <- const_df()
    suffix <- if (input$textPosition10 == "up") "text10up" else "text10down"
    varName <- if (input$level10 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    
    df <- df[df$Var != varName, ]
    new_row <- data.frame(Var = varName, Select = input$text10, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
    
    updateTextAreaInput(session, "text10", value = "")
  })
  
  observeEvent(input$deleteText10, {
    df <- const_df()
    suffix <- if (input$textPosition10 == "up") "text10up" else "text10down"
    varName <- if (input$level10 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    
    df <- df[df$Var != varName, ]
    const_df(df)
    
    updateTextAreaInput(session, "text10", value = "")
  })
  
  observe({
    df <- const_df()
    suffix <- if (input$textPosition10 == "up") "text10up" else "text10down"
    varName <- if (input$level10 == "Protein") {
      paste0(suffix, "Protein")
    } else {
      paste0(suffix, "Phosphosite")
    }
    
    value <- df$Select[df$Var == varName]
    if (length(value) > 0) {
      updateTextAreaInput(session, "text10", value = value)
    } else {
      updateTextAreaInput(session, "text10", value = "")
    }
  })
  
  output$table10 <- DT::renderDataTable({
    req(data10())
    DT::datatable(data10(), options = list(pageLength = 20), callback = JS("
      $.fn.dataTable.ext.errMode = 'none';
      $(document).on('error.dt', function(e, settings, techNote, message) {
        console.warn('Suppressed DataTables warning:', message);
        e.preventDefault();
      });
    "))
  })
  
  output$download10<- downloadHandler(
    filename = function() {
      paste("volcanodata_", Sys.Date(), "_", input$condition1_10 ,"_", input$condition2_10 ,".csv", sep = "")
    },
    content = function(file) {
      write.csv(data10(), file, row.names = FALSE)
    }
  )
  
  output$download10a <- downloadHandler(
    filename = function() {
      paste0("volcano_all_comparisons_", format(Sys.Date(), "%Y%m%d"), ".xlsx")
    },
    content = function(file) {
      file_path <- all_volcano_data(
        data = volcano_data(),       
        meta = meta(),                 
        in_pval = input$in_pval10,     
        in_log2fc = input$in_log2fc10,
        workflow = input$level10,      
        paired = input$paired10      
      )
      
      file.copy(file_path, file, overwrite = TRUE)
    }
  )
  
#### Tab11 - GSEA ####
  different_genes_df <- reactive({
    req(volcano_data(), meta(), input$condition1_10, input$condition2_10)
    different_genes(volcano_data(), meta(), input$condition1_10, input$condition2_10, in_pval = input$in_pval10, in_log2fc = input$in_log2fc10)
  })
  
  output$UpregEnrichmentPlot <- renderPlot({
    req(different_genes_df())
    enrichment_analysis(different_genes_df()$Upregulated, top_n = input$top_n11, min_num = input$filter11min, max_num = input$filter11max)
  })
  
  output$DownregEnrichmentPlot <- renderPlot({
    req(different_genes_df())
    enrichment_analysis(different_genes_df()$Downregulated, top_n = input$top_n11, min_num = input$filter11min, max_num = input$filter11max)
  })
  
  output$downloadUpPlot <- downloadHandler(
    filename = function() {
      paste("enrichup_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth11
      height <- input$plotHeight11
      dpi <- input$plotDPI11
      
      ggsave(file, plot = enrichment_analysis(different_genes_df()$Upregulated, top_n = input$top_n11, min_num = input$filter11min, max_num = input$filter11max), device = "png", 
             width = width, height = height, units = "cm", dpi = dpi)
    }
  )
  
  output$downloadDownPlot <- downloadHandler(
    filename = function() {
      paste("enrichdown_plot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth11
      height <- input$plotHeight11
      dpi <- input$plotDPI11
      
      ggsave(file, plot = enrichment_analysis(different_genes_df()$Downregulated, top_n = input$top_n11, min_num = input$filter11min, max_num = input$filter11max), device = "png", 
             width = width, height = height, units = "cm", dpi = dpi)
    }
  )
  
  output$table11up <- DT::renderDataTable({
    req(different_genes_df())
    DT::datatable(different_genes_df()$Upregulated, options = list(pageLength = 5))
  })
  
  output$table11down <- DT::renderDataTable({
    req(different_genes_df())
    DT::datatable(different_genes_df()$Downregulated, options = list(pageLength = 5))
  })
  
#### Tab10.5 - Simulation ####  
  observe({
    if (input$level10.5 == "Protein") {
      req(transformed_data())
      updateSelectInput(session, "condition1_10.5", choices = unique(meta()$condition))
      updateSelectInput(session, "condition2_10.5", choices = unique(meta()$condition))
    } else if (input$level10.5 == "Phosphosite") {
      req(transformed_data3())
      updateSelectInput(session, "condition1_10.5", choices = unique(meta2()$condition))
      updateSelectInput(session, "condition2_10.5", choices = unique(meta2()$condition))
    }
  })
  
  output$VolcPlotSim <- renderPlotly({
    req(input$condition1_10.5, input$condition2_10.5)
    
    if (input$level10.5 == "Protein") {
      req(volcano_data(), meta())
      volcano_plot_sim(
        data       = volcano_data(),
        meta       = meta(),
        condition1 = input$condition1_10.5,
        condition2 = input$condition2_10.5,
        in_pval    = input$in_pval10.5,
        in_log2fc  = input$in_log2fc10.5,
        mod_var    = input$mod_var10.5,
        mod_n      = ifelse(input$mod_n10.5 == 1, 0, input$mod_n10.5),  # 0 means use real n
        workflow   = "Protein"
      )
    } else if (input$level10.5 == "Phosphosite") {
      req(volcano_data3(), meta2())
      volcano_plot_sim(
        data       = volcano_data3(),
        meta       = meta2(),
        condition1 = input$condition1_10.5,
        condition2 = input$condition2_10.5,
        in_pval    = input$in_pval10.5,
        in_log2fc  = input$in_log2fc10.5,
        mod_var    = input$mod_var10.5,
        mod_n      = ifelse(input$mod_n10.5 == 1, 0, input$mod_n10.5),
        workflow   = "Phosphosite"
      )
    }
  })
  
  
#### Tab12 - Correlation Plot ####
  id12 <- reactiveVal(TRUE)
  MeCorr <- reactiveVal(FALSE)
  
  observeEvent(input$toggle_id12, {
    id12(!id12())
  })
  
  observeEvent(input$Change12, {
    MeCorr(!MeCorr())
  })
  
  output$CorrPlot <- renderPlot({
    if (input$level12 == "Protein") {
      req(data(), meta())
      corr_plot(data(), meta(), MeCorr(), id12())
    } else if (input$level12 == "Phosphosite") {
      req(data3(), meta2())
      corr_plot(data3(), meta2(), MeCorr(), id12())
    }
  })
  
  observe({
    df <- const_df()
    text_up <- ""
    text_down <- ""
    
    if (input$level12 == "Protein") {
      if ("text12upProtein" %in% df$Var) {
        text_up <- df$Select[df$Var == "text12upProtein"]
      }
      if ("text12downProtein" %in% df$Var) {
        text_down <- df$Select[df$Var == "text12downProtein"]
      }
    } else if (input$level12 == "Phosphosite") {
      if ("text12upPhosphosite" %in% df$Var) {
        text_up <- df$Select[df$Var == "text12upPhosphosite"]
      }
      if ("text12downPhosphosite" %in% df$Var) {
        text_down <- df$Select[df$Var == "text12downPhosphosite"]
      }
    }
    
    if (input$textPosition12 == "up") {
      updateTextAreaInput(session, "text12", value = text_up)
    } else {
      updateTextAreaInput(session, "text12", value = text_down)
    }
  })
  
  observeEvent(input$addText12, {
    df <- const_df()
    
    if (input$level12 == "Protein") {
      if (input$textPosition12 == "up") {
        df <- df[df$Var != "text12upProtein", ]
        new_row <- data.frame(Var = "text12upProtein", Select = input$text12, stringsAsFactors = FALSE)
      } else {
        df <- df[df$Var != "text12downProtein", ]
        new_row <- data.frame(Var = "text12downProtein", Select = input$text12, stringsAsFactors = FALSE)
      }
    } else if (input$level12 == "Phosphosite") {
      if (input$textPosition12 == "up") {
        df <- df[df$Var != "text12upPhosphosite", ]
        new_row <- data.frame(Var = "text12upPhosphosite", Select = input$text12, stringsAsFactors = FALSE)
      } else {
        df <- df[df$Var != "text12downPhosphosite", ]
        new_row <- data.frame(Var = "text12downPhosphosite", Select = input$text12, stringsAsFactors = FALSE)
      }
    }
    
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text12", value = "")
  })
  
  observeEvent(input$deleteText12, {
    df <- const_df()
    
    if (input$level12 == "Protein") {
      if (input$textPosition12 == "up") {
        df <- df[df$Var != "text12upProtein", ]
      } else {
        df <- df[df$Var != "text12downProtein", ]
      }
    } else if (input$level12 == "Phosphosite") {
      if (input$textPosition12 == "up") {
        df <- df[df$Var != "text12upPhosphosite", ]
      } else {
        df <- df[df$Var != "text12downPhosphosite", ]
      }
    }
    
    const_df(df)
    updateTextAreaInput(session, "text12", value = "")
  })
  
    
#### Tab13 - Heatmap ####
  id13 <- reactiveVal(TRUE)
  show13 <- reactiveVal(FALSE)
  
  observeEvent(input$toggle_id13, {
    id13(!id13())
  })
  
  observeEvent(input$show13, {
    show13(!show13())
  })
  
  observeEvent(input$fix13, {
    clear_all_plots()
  })
  
  output$heatPlot <- renderPlot({
    if (input$level13=="Protein"){
      req(data(), meta())
      if (show13()==TRUE){
        heatmap_plot(data(), meta(), id13())
      } else if (show13()==FALSE){
        heatmap_plot_nmv(transformed_data(), meta(), id13())
        }
    } else if (input$level13=="Phosphosite"){
      req(data3(), meta2())
      if (show13()==TRUE){
        heatmap_plot(data3(), meta2(), id13())
      } else if (show13()==FALSE){
        heatmap_plot_nmv(transformed_data3(), meta2(), id13())
      }
    }
  })
  
  output$downloadheatPlot <- downloadHandler(
    filename = function() {
      paste("heatPlot", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      width <- input$plotWidth13
      height <- input$plotHeight13 / 2.54 
      dpi <- input$plotDPI13 / 2.54 
      
      png(file, width = width, height = height, units = "in", res = dpi)
      heatmap_plot(data(), meta(), id13())
      dev.off()
    }
  )
  
  observe({
    df <- const_df()
    text_up <- ""
    text_down <- ""
    
    if (input$level13 == "Protein") {
      if ("text13upProtein" %in% df$Var) {
        text_up <- df$Select[df$Var == "text13upProtein"]
      }
      if ("text13downProtein" %in% df$Var) {
        text_down <- df$Select[df$Var == "text13downProtein"]
      }
    } else if (input$level13 == "Phosphosite") {
      if ("text13upPhosphosite" %in% df$Var) {
        text_up <- df$Select[df$Var == "text13upPhosphosite"]
      }
      if ("text13downPhosphosite" %in% df$Var) {
        text_down <- df$Select[df$Var == "text13downPhosphosite"]
      }
    }
    
    if (input$textPosition13 == "up") {
      updateTextAreaInput(session, "text13", value = text_up)
    } else {
      updateTextAreaInput(session, "text13", value = text_down)
    }
  })
  
  observeEvent(input$addText13, {
    df <- const_df()
    
    if (input$level13 == "Protein") {
      if (input$textPosition13 == "up") {
        df <- df[df$Var != "text13upProtein", ]
        new_row <- data.frame(Var = "text13upProtein", Select = input$text13, stringsAsFactors = FALSE)
      } else {
        df <- df[df$Var != "text13downProtein", ]
        new_row <- data.frame(Var = "text13downProtein", Select = input$text13, stringsAsFactors = FALSE)
      }
    } else if (input$level13 == "Phosphosite") {
      if (input$textPosition13 == "up") {
        df <- df[df$Var != "text13upPhosphosite", ]
        new_row <- data.frame(Var = "text13upPhosphosite", Select = input$text13, stringsAsFactors = FALSE)
      } else {
        df <- df[df$Var != "text13downPhosphosite", ]
        new_row <- data.frame(Var = "text13downPhosphosite", Select = input$text13, stringsAsFactors = FALSE)
      }
    }
    
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text13", value = "")
  })
  
  observeEvent(input$deleteText13, {
    df <- const_df()
    
    if (input$level13 == "Protein") {
      if (input$textPosition13 == "up") {
        df <- df[df$Var != "text13upProtein", ]
      } else {
        df <- df[df$Var != "text13downProtein", ]
      }
    } else if (input$level13 == "Phosphosite") {
      if (input$textPosition13 == "up") {
        df <- df[df$Var != "text13upPhosphosite", ]
      } else {
        df <- df[df$Var != "text13downPhosphosite", ]
      }
    }
    
    const_df(df)
    updateTextAreaInput(session, "text13", value = "")
  })
  
#### Tab14 - RT Plot ####
  line14 <- reactiveVal(FALSE)
  header14 <- reactiveVal(TRUE)
  
  observeEvent(input$line14, {
    line14(!line14())
  })
  
  observeEvent(input$toggle_header14, {
    header14(!header14())
  })

  output$RTPlot <- renderPlot({
    req(data2(), meta())
    RTvspredRT_plot(
      data2 = data2(),
      meta = meta(),
      method = input$style14,
      add_line = line14(),
      bin = input$bins14,
      header = header14()
    )
  })
  
  output$downloadRTPlot <- downloadHandler(
    filename = function() {
      paste0("RT_plot_", Sys.Date(), ".", input$plotFormat14)
    },
    content = function(file) {
      req(data2(), meta())
      
      p <- RTvspredRT_plot(
        data2 = data2(),
        meta = meta(),
        method = input$style14,
        add_line = line14(),
        bin = input$bins14,
        header = header14()
      )
      
      width <- input$plotWidth14
      height <- input$plotHeight14
      dpi <- input$plotDPI14
      ext <- input$plotFormat14
      
      if (ext %in% c("png", "jpg")) {
        ggsave(file, plot = p, device = ext, width = width, height = height, units = "cm", dpi = dpi)
      } else if (ext == "pdf") {
        pdf(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      } else if (ext == "svg") {
        svglite::svglite(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      }
    }
  )
  
  observe({
    df <- const_df()
    text_up <- ""
    text_down <- ""
    
    if ("text14up" %in% df$Var) {
      text_up <- df$Select[df$Var == "text14up"]
    }
    if ("text14down" %in% df$Var) {
      text_down <- df$Select[df$Var == "text14down"]
    }
    
    if (input$textPosition14 == "up") {
      updateTextAreaInput(session, "text14", value = text_up)
    } else {
      updateTextAreaInput(session, "text14", value = text_down)
    }
  })
  
  observeEvent(input$addText14, {
    df <- const_df()
    
    if (input$textPosition14 == "up") {
      df <- df[df$Var != "text14up", ]
      new_row <- data.frame(Var = "text14up", Select = input$text14, stringsAsFactors = FALSE)
    } else {
      df <- df[df$Var != "text14down", ]
      new_row <- data.frame(Var = "text14down", Select = input$text14, stringsAsFactors = FALSE)
    }
    
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text14", value = "")
  })
  
  observeEvent(input$deleteText14, {
    df <- const_df()
    
    if (input$textPosition14 == "up") {
      df <- df[df$Var != "text14up", ]
    } else {
      df <- df[df$Var != "text14down", ]
    }
    
    const_df(df)
    updateTextAreaInput(session, "text14", value = "")
  })
  
#### Tab15 - Modification Plot ####
  id15 <- reactiveVal(TRUE)
  header15 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_id15, {
    id15(!id15())
  })
  
  observeEvent(input$toggle_header15, {
    header15(!header15())
  })
  
  output$ModPlot <- renderPlot({
    req(data2(), meta())
    modification_plot(data2 = data2(), meta = meta(), id = id15(), header = header15())
  })
  
  output$downloadModPlot <- downloadHandler(
    filename = function() {
      ext <- input$plotFormat15
      paste("modification_plot", Sys.Date(), ".", ext, sep = "")
    },
    content = function(file) {
      req(data2(), meta())
      p <- modification_plot(data2 = data2(), meta = meta(), id = id15(), header = header15())
      width <- input$plotWidth15
      height <- input$plotHeight15
      dpi <- input$plotDPI15
      ext <- input$plotFormat15
      
      if (ext %in% c("png", "jpg")) {
        ggsave(file, plot = p, device = ext, width = width, height = height, units = "cm", dpi = dpi)
      } else if (ext == "pdf") {
        pdf(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      } else if (ext == "svg") {
        svglite::svglite(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      }
    }
  )
  
  observe({
    df <- const_df()
    text_up <- ""
    text_down <- ""
    
    if ("text15up" %in% df$Var) {
      text_up <- df$Select[df$Var == "text15up"]
    }
    if ("text15down" %in% df$Var) {
      text_down <- df$Select[df$Var == "text15down"]
    }
    
    if (input$textPosition15 == "up") {
      updateTextAreaInput(session, "text15", value = text_up)
    } else {
      updateTextAreaInput(session, "text15", value = text_down)
    }
  })
  
  observeEvent(input$addText15, {
    df <- const_df()
    
    if (input$textPosition15 == "up") {
      df <- df[df$Var != "text15up", ]
      new_row <- data.frame(Var = "text15up", Select = input$text15, stringsAsFactors = FALSE)
    } else {
      df <- df[df$Var != "text15down", ]
      new_row <- data.frame(Var = "text15down", Select = input$text15, stringsAsFactors = FALSE)
    }
    
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text15", value = "")
  })
  
  observeEvent(input$deleteText15, {
    df <- const_df()
    
    if (input$textPosition15 == "up") {
      df <- df[df$Var != "text15up", ]
    } else {
      df <- df[df$Var != "text15down", ]
    }
    
    const_df(df)
    updateTextAreaInput(session, "text15", value = "")
  })
  
#### Tab16 - Missed Cleavage Plot ####
  id16 <- reactiveVal(TRUE)
  header16 <- reactiveVal(TRUE)
  text16 <- reactiveVal(TRUE)
  
  observeEvent(input$toggle_id16, {
    id16(!id16())
  })
  
  observeEvent(input$toggle_header16, {
    header16(!header16())
  })
  
  observeEvent(input$toggle_text16, {
    text16(!text16())
  })
  
  output$MCPlot <- renderPlot({
    req(data2(), meta())
    missed_cl_plot(data2 = data2(), meta = meta(), id = id16(), text = text16(), text_size = input$text_size16, header = header16())
  })
  
  output$downloadMCPlot <- downloadHandler(
    filename = function() {
      ext <- input$plotFormat16
      paste("missed_cleavage_plot", Sys.Date(), ".", ext, sep = "")
    },
    content = function(file) {
      req(data2(), meta())
      p <- missed_cl_plot(data2 = data2(), meta = meta(), id = id16(), text = text16(), text_size = input$text_size16, header = header16())
      width <- input$plotWidth16
      height <- input$plotHeight16
      dpi <- input$plotDPI16
      ext <- input$plotFormat16
      
      if (ext %in% c("png", "jpg")) {
        ggsave(file, plot = p, device = ext, width = width, height = height, units = "cm", dpi = dpi)
      } else if (ext == "pdf") {
        pdf(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      } else if (ext == "svg") {
        svglite::svglite(file, width = width / 2.54, height = height / 2.54)
        print(p)
        dev.off()
      }
    }
  )
  
  observe({
    df <- const_df()
    text_up <- ""
    text_down <- ""
    
    if ("text16up" %in% df$Var) {
      text_up <- df$Select[df$Var == "text16up"]
    }
    if ("text16down" %in% df$Var) {
      text_down <- df$Select[df$Var == "text16down"]
    }
    
    if (input$textPosition16 == "up") {
      updateTextAreaInput(session, "text16", value = text_up)
    } else {
      updateTextAreaInput(session, "text16", value = text_down)
    }
  })
  
  observeEvent(input$addText16, {
    df <- const_df()
    
    if (input$textPosition16 == "up") {
      df <- df[df$Var != "text16up", ]
      new_row <- data.frame(Var = "text16up", Select = input$text16, stringsAsFactors = FALSE)
    } else {
      df <- df[df$Var != "text16down", ]
      new_row <- data.frame(Var = "text16down", Select = input$text16, stringsAsFactors = FALSE)
    }
    
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text16", value = "")
  })
  
  observeEvent(input$deleteText16, {
    df <- const_df()
    
    if (input$textPosition16 == "up") {
      df <- df[df$Var != "text16up", ]
    } else {
      df <- df[df$Var != "text16down", ]
    }
    
    const_df(df)
    updateTextAreaInput(session, "text16", value = "")
  })
  
#### Tab17 - Boxplot, single protein ####
  observe({
    if (input$level17 == "Protein"){
      req(transformed_data())
      updateSelectizeInput(session, "protein17", choices = unique(data()$ProteinNames))
      updateSelectizeInput(session, "conditions17", choices = unique(meta()$condition))
    } else if (input$level17 == "Phosphosite"){
      req(transformed_data())
      updateSelectizeInput(session, "protein17", choices = unique(data3()$PTM_Collapse_key))
      updateSelectizeInput(session, "conditions17", choices = unique(meta2()$condition))
    }
  })
  
  output$ProtBoxPlot <- renderPlot({
    if (input$level17 == "Protein"){
      req(transformed_data(), meta())
      validate(
        need(length(input$conditions17) >= 2, "Please select at least two conditions")
      )
      plot <- tryCatch({
        compare_prot_box(transformed_data(), meta(), input$conditions17, input$protein17)
      }, error = function(e) {
        stop(safeError(e))
      })
      print(plot)
    } else if (input$level17=="Phosphosite"){
      req(transformed_data3(), meta2())
      validate(
        need(length(input$conditions17) >= 2, "Please select at least two conditions")
      )
      plot <- tryCatch({
        compare_prot_box(transformed_data3(), meta2(), input$conditions17, input$protein17, workflow = "Phosphosite")
      }, error = function(e) {
        stop(safeError(e))
      })
      print(plot)
    }
  })
  
#### Tab18 - Lineplot, single protein ####
  observe({
    if (input$level18 == "Protein"){
      req(transformed_data())
      updateSelectizeInput(session, "protein18", choices = unique(data()$ProteinNames))
      updateSelectizeInput(session, "conditions18", choices = unique(meta()$condition))
    } else if (input$level18 == "Phosphosite"){
      req(transformed_data())
      updateSelectizeInput(session, "protein18", choices = unique(data3()$PTM_Collapse_key))
      updateSelectizeInput(session, "conditions18", choices = unique(meta2()$condition))
    }
  })
  
  output$ProtLinePlot <- renderPlot({
    if (input$level18=="Protein"){
      req(transformed_data(), meta())
      plot <- tryCatch({
        compare_prot_line(transformed_data(), meta(), input$conditions18, input$protein18)
      }, error = function(e) {
        stop(safeError(e))
      })
      print(plot)
    } else if (input$level18=="Phosphosite"){
      req(transformed_data3(), meta2())
      plot <- tryCatch({
        compare_prot_line(transformed_data3(), meta2(), input$conditions18, input$protein18, workflow = "Phosphosite")
      }, error = function(e) {
        stop(safeError(e))
      })
      print(plot)
    }
  })
  
#### Tab19 - Phossite Plot ####
  output$PhossitePlot <- renderPlot({
    req(data3())
    simple_phos_site_plot(data3(), filter=input$cutoff19)
  })
  
  observe({
    df <- const_df()
    text_up <- ""
    text_down <- ""
    
    if ("text19up" %in% df$Var) {
      text_up <- df$Select[df$Var == "text19up"]
    }
    if ("text19down" %in% df$Var) {
      text_down <- df$Select[df$Var == "text19down"]
    }
    
    if (input$textPosition19 == "up") {
      updateTextAreaInput(session, "text19", value = text_up)
    } else {
      updateTextAreaInput(session, "text19", value = text_down)
    }
  })
  
  observeEvent(input$addText19, {
    df <- const_df()
    
    if (input$textPosition19 == "up") {
      df <- df[df$Var != "text19up", ]
      new_row <- data.frame(Var = "text19up", Select = input$text19, stringsAsFactors = FALSE)
    } else {
      df <- df[df$Var != "text19down", ]
      new_row <- data.frame(Var = "text19down", Select = input$text19, stringsAsFactors = FALSE)
    }
    
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text19", value = "")
  })
  
  observeEvent(input$deleteText19, {
    df <- const_df()
    
    if (input$textPosition19 == "up") {
      df <- df[df$Var != "text19up", ]
    } else {
      df <- df[df$Var != "text19down", ]
    }
    
    const_df(df)
    updateTextAreaInput(session, "text19", value = "")
  })
  
#### Tab20.1 - KSEA ####
  observe({
    req(transformed_data3())
    updateSelectInput(session, "condition1_20", 
                      choices = c(unique(meta2()$condition)))
  })
  
  observe({
    req(transformed_data3())
    updateSelectInput(session, "condition2_20", 
                      choices = c(unique(meta2()$condition)))
  })
  
  output$downKSEA <- downloadHandler(
    filename = function() {
      paste("KSEA_data_", input$condition1_20, "_", input$condition2_20, ".tsv", sep = "")  
    },
    content = function(file) {
      ksea_result <- prepare_KSEA(
        data = volcano_data3(),
        meta = meta2(),
        condition1 = input$condition1_20,
        condition2 = input$condition2_20,
        pvalue = input$include_pval_ksea
      )
      
      write.table(
        ksea_result,
        file,
        sep = "\t",
        row.names = FALSE,
        quote = FALSE
      )
    }
  )
  
  KSEA_data <- reactive({
    req(input$KSEA_data)
    file <- input$KSEA_data$datapath
    read.delim(file, sep = "\t")
  })
  
  output$KinaseVolcPlot <- renderPlotly({
    req(KSEA_data(), input$in_pval20)
    kinase_volcano(KSEA_data(), input$in_pval20)
  })
  
  output$downTREE <- downloadHandler(
    filename = function() {
      paste("TREE_data_", input$condition1_20, "_", input$condition2_20, ".zip", sep = "")  
    },
    content = function(file) {
      temp_dir <- tempdir()
      
      file1 <- file.path(temp_dir, "enrichment_scores.txt")
      file2 <- file.path(temp_dir, "pvals.txt")
      
      write.table(prepare_tree_data_es(KSEA_data()), 
                  file1, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
      
      write.table(prepare_tree_data_pv(KSEA_data()), 
                  file2, sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)
      
      zip::zipr(file, files = c(file1, file2))
    }
  )
  
  output$treePlot <- renderUI({
    req(input$TreeSVG)
    svg_content <- readLines(input$TreeSVG$datapath)
    HTML(paste(svg_content, collapse = "\n"))
  })
  
#### Tab20.2 - KSEA (Kinact)####
  observe({
    req(transformed_data3())
    updateSelectInput(session, "condition1_202", 
                      choices = c(unique(meta2()$condition)))
  })
  
  observe({
    req(transformed_data3())
    updateSelectInput(session, "condition2_202", 
                      choices = c(unique(meta2()$condition)))
  })
  
  data202 <- reactive({
    req(volcano_data3(), meta2(), input$condition1_202, input$condition2_202)
    volc_data <- volcano_data_f(volcano_data3(), meta2(), input$condition1_202, input$condition2_202, workflow="Phosphosite")
    rownames(volc_data) <- NULL
    org_data = volcano_data3()
    data = match_volc_and_org_data(volc_data = volc_data, org_data = org_data)
    return(data)
  })
  
  output$KinactPlot <- renderPlotly({
      req(data202())
      kinact_kinase_activity(data202(), top_n = input$top_n202, NetworKIN = input$NetworKIN202, NetworKIN.cutoff = input$NetworKIN.cutoff202, m.cutoff = input$m.cutoff202)
    })
  
  output$KinactPlotSites <- renderPlotly({
    req(data202(), input$Kinase_202)
    downstream_phossite_volc(data202(), input$Kinase_202, NetworKIN = input$NetworKIN202)
  })
  
#### Tab21 - Phossite Coverage Plot ####
  output$PhossiteCoveragePlot <- renderPlot({
    req(data3(), meta2())
    phossite_coverage_plot(data3(), meta2())
  })
  
  observe({
    df <- const_df()
    text_up <- ""
    text_down <- ""
    
    if ("text21up" %in% df$Var) {
      text_up <- df$Select[df$Var == "text21up"]
    }
    if ("text21down" %in% df$Var) {
      text_down <- df$Select[df$Var == "text21down"]
    }
    
    if (input$textPosition21 == "up") {
      updateTextAreaInput(session, "text21", value = text_up)
    } else {
      updateTextAreaInput(session, "text21", value = text_down)
    }
  })
  
  observeEvent(input$addText21, {
    df <- const_df()
    
    if (input$textPosition21 == "up") {
      df <- df[df$Var != "text21up", ]
      new_row <- data.frame(Var = "text21up", Select = input$text21, stringsAsFactors = FALSE)
    } else {
      df <- df[df$Var != "text21down", ]
      new_row <- data.frame(Var = "text21down", Select = input$text21, stringsAsFactors = FALSE)
    }
    
    const_df(rbind(df, new_row))
    updateTextAreaInput(session, "text21", value = "")
  })
  
  observeEvent(input$deleteText21, {
    df <- const_df()
    
    if (input$textPosition21 == "up") {
      df <- df[df$Var != "text21up", ]
    } else {
      df <- df[df$Var != "text21down", ]
    }
    
    const_df(df)
    updateTextAreaInput(session, "text21", value = "")
  })
  
#### Tab22 - Protein centric ####
  db_fasta <- reactiveVal(NULL)
  
  observeEvent(input$db_fasta22, {
    #fasta_path <- file.path(getwd(), "www", "db", "UP000005640_9606.fasta") 
    #fasta_path <- file.path(getwd(), "www", "db", "UP000005640_9606.fasta") 
    fasta_path <- "C:/Users/Jonas Marx/Desktop/Proteomics CopilotR/www/db/UP000000589_10090.fasta"
    db_fasta_stat <- readFASTA(fasta_path)
    db_fasta_stat = as.data.frame(db_fasta_stat)
    db_fasta_stat = as.data.frame(t(db_fasta_stat))
    db_fasta_stat$name <- rownames(db_fasta_stat)
    db_fasta_stat$name <- sapply(strsplit(as.character(db_fasta_stat$name), "\\."), tail, 1)
    db_fasta(db_fasta_stat)
    
    showNotification("Database loaded successfully!", type = "message", duration = 4)
  })
  
  observe({
    req(transformed_data())
    updateSelectizeInput(session, "protein20", choices = unique(data()$ProteinNames))
  })
  
  output$seq_allign <- renderPrint({
    split_string <- strsplit(input$protein20, ";")[[1]]
    vis_coverage(data2(), split_string[1], chunk_size = input$chunk20, db = db_fasta())
  })
  
  output$plot_3d22 <- renderR3dmol({
    protein <- input$protein20
    plot <- model_3d(data2(), protein, db = db_fasta())
    return(plot)
  })
  
  
#### Tabn1 - Impute data ####
  output$Impute1 <- renderPlot({
    if (input$leveln1 == "Protein") {
      req(transformed_data(), meta())
      impute_values(transformed_data(), meta(), ret=1, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    } else if (input$leveln1 == "Phosphosite"){
      req(transformed_data3(), meta2())
      impute_values(transformed_data3(), meta2(), ret=1, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    }
  })
  
  output$Impute2 <- renderPlot({
    if (input$leveln1 == "Protein") {
      req(transformed_data(), meta())
      impute_values(transformed_data(), meta(), ret=2, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    } else if (input$leveln1 == "Phosphosite"){
      req(transformed_data3(), meta2())
      impute_values(transformed_data3(), meta2(), ret=2, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    }
  })
  
  output$Impute3 <- renderPlot({
    if (input$leveln1 == "Protein") {
      req(transformed_data(), meta())
      impute_values(transformed_data(), meta(), ret=3, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    } else if (input$leveln1 == "Phosphosite"){
      req(transformed_data3(), meta2())
      impute_values(transformed_data3(), meta2(), ret=3, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1)
    }
  })
  
  was_imputed <- reactiveVal(FALSE)
  was_imputed3 <- reactiveVal(FALSE)
  imputed_data <- reactiveVal(NULL)
  imputed_data3 <- reactiveVal(NULL)
  
  observeEvent(input$ImputeEVE, {
    if (input$leveln1 == "Protein") {
      req(pre_volcano_data(), meta())
      imputed_data(impute_values(pre_volcano_data(), meta(), ret=0, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1))
      volcano_data <- reactive({
        imputed_data(TRUE)
      })
      was_imputed(TRUE)
    } else if (input$leveln1 == "Phosphosite"){
      req(pre_volcano_data3(), meta2())
      imputed_data3(impute_values(pre_volcano_data3(), meta2(), ret=0, q=input$qn1, adj_std = input$adj_stdn1, seed=input$seedn1))
      volcano_data3 <- reactive({
        imputed_data3()
      })
      was_imputed3(TRUE)
    }
  })
  
  observeEvent(input$ImputeEVE, {
    if (input$leveln1 == "Protein") {
      req(volcano_data(), meta(), input$qn1, input$adj_stdn1, input$seedn1)
      df <- const_df()
      df <- df[!df$Var %in% c("Q-shift", "Adj-std", "Seed"), ]
      new_rows <- data.frame(
        Var = c("Q-shift", "Adj-std", "Seed"),
        Select = c(input$qn1, input$adj_stdn1, input$seedn1),
        stringsAsFactors = FALSE
      )
      const_df(rbind(df, new_rows))
    } else if (input$leveln1 == "Phosphosite"){
      req(volcano_data3(), meta2(), input$qn1, input$adj_stdn1, input$seedn1)
      df <- const_df()
      df <- df[!df$Var %in% c("Q-shift2", "Adj-std2", "Seed2"), ]
      new_rows <- data.frame(
        Var = c("Q-shift2", "Adj-std2", "Seed2"),
        Select = c(input$qn1, input$adj_stdn1, input$seedn1),
        stringsAsFactors = FALSE
      )
      const_df(rbind(df, new_rows))
    }
  })
  
  output$imputed_data_tab <- DT::renderDataTable({
    if (input$leveln1 == "Protein") {
      req(imputed_data())
      DT::datatable(imputed_data(), options = list(pageLength = 10))
    } else if (input$leveln1 == "Phosphosite"){
      req(imputed_data3())
      DT::datatable(imputed_data3(), options = list(pageLength = 10))
    }
  })
  
  output$imputed_data_down <- downloadHandler(
      filename = function() {
        if (input$leveln1 == "Protein") {
          paste("imputed_data-", Sys.Date(), ".csv", sep = "")
        } else if (input$leveln1 == "Phosphosite") {
          paste("imputed_data_phospho-", Sys.Date(), ".csv", sep = "")
        }
      },
      content = function(file) {
        if (input$leveln1 == "Protein") {
          req(imputed_data())
          write.csv(imputed_data(), file, row.names = FALSE)
        } else if (input$leveln1 == "Phosphosite") {
          req(imputed_data3())
          write.csv(imputed_data3(), file, row.names = FALSE)
        }
      }
  )
  
#### TabSummary - Summary ####
  output$transformed_data_prot <- DT::renderDataTable({
    req(transformed_data())
    DT::datatable(transformed_data(), options = list(pageLength = 5))
  })
  
  output$meta_object <- DT::renderDataTable({
    req(meta())
    DT::datatable(meta(), options = list(pageLength = 5))
  })
  
  output$collapse_data_sum <- DT::renderDataTable({
    req(data3())
    DT::datatable(data3(), options = list(pageLength = 5))
  })
  
  output$meta_object2 <- DT::renderDataTable({
    req(meta2())
    DT::datatable(meta2(), options = list(pageLength = 5))
  })
  
  output$download_meta <- downloadHandler(
    filename = function() {
      paste("metadata-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(meta(), file, row.names = FALSE)
    }
  )
  
  output$download_meta2 <- downloadHandler(
    filename = function() {
      paste("metadata_phos-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(meta2(), file, row.names = FALSE)
    }
  )
  
  output$download_collapse_data_sum <- downloadHandler(
    filename = function() {
      paste("collapse_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      write.csv(data3(), file, row.names = FALSE)
    }
  )
  
#### LogLogic ####
  rv <- reactiveValues(data = NULL)
  
  observe({
    var_list <- list(
      #Version
      "Version",
      
      #Color
      "Color",
      
      #Meta
      "Transformed", "Transformed3",
      
      #Imputation
      "Imputed", "Imputed3",
      
      #Coverage Plot
      "CoveragePlotID", "CoveragePlotHeader", "CoveragePlotLegend",
      
      #Missing Value Plot
      "MissingValuePlotHeader", "MissingValuePlotBin", "MissingValueText", "MissingValueTextSize",
      
      #Histogram Intensity
      "HistogramIntHeader", "HistogramIntLegend",
      
      #Boxplot Intensity
      "BoxplotIntID", "BoxplotIntOut", "BoxplotIntMean", "BoxplotIntHeader", "BoxplotIntLegend",
      
      #COV Plot
      "CovPlotOut", "CovPlotHeader", "CovPlotLegend",
      
      #PCA Plot
      "PCAHeader", "PCALegend",
      
      #CorrPlot
      "CorrPlotDisplay", "CorrPlotID",
      
      #Heatmap
      "HeatmapID",
      
      #RT Plot
      "RTPlotSytle", "RTline", "HexbinsRT", "RTHeader",
      
      #Modifications Plot
      "ModPlotID", "ModPlotHeader",
      
      #Missied Cleavage Plot
      "MissedCleavID", "MissedCleavText", "MissedCleavTextSize", "MissedCleavHeader")
    
    select_list <- list(
      #Version
      "V0.1.1",
      
      #Color
      input$color_palette,
      
      #Meta
      was_transformed(), was_transformed3(),
      
      #Imputation
      was_imputed(), was_imputed3(),
      
      #Coverage Plot
      id3(), header3(), legend3(),
      
      #Missing Value Plot
      header4(), input$missValBin4, text4(), input$text_size4,
      
      #Histogram Intensity
      header5(), legend5(),
      
      #Boxplot Intensity
      id6(), outliersBX(), meanBX(), header6(), legend6(),
      
      #COV Plot
      outliersCOV(), header7(), legend7(),
      
      #PCA Plot
      header8(), legend8(),
      
      #CorrPlot
      MeCorr(), id12(), 
      
      #Heatmap
      id13(), 
      
      #RT Plot
      input$style14, line14(), input$bins14, header14(),
      
      #Modifications Plot
      id15(), header15(),
      
      #Missed Cleavage Plot
      id16(), text16(), input$text_size16, header16())
    
    df <- do.call(rbind, Map(data.frame, Var = var_list, Select = select_list, stringsAsFactors = FALSE))
    df$Select[df$Select == 1] <- TRUE
    df$Select[df$Select == 0] <- FALSE
    rv$data <- df
  })
  
  const_df <- reactiveVal(data.frame(Var = character(), Select = character(), stringsAsFactors = FALSE))
  
  observeEvent(input$file, {
    req(input$file)
    df <- const_df()
    file_name <- input$file$name
    df <- df[df$Var != "Data1", ]
    new_row <- data.frame(Var = "Data1", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$file2, {
    req(input$file2)
    df <- const_df()
    file_name <- input$file2$name
    df <- df[df$Var != "Data2", ]
    new_row <- data.frame(Var = "Data2", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$file3, {
    req(input$file3)
    df <- const_df()
    file_name <- input$file3$name
    df <- df[df$Var != "Data3", ]
    new_row <- data.frame(Var = "Data3", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$KSEA_data, {
    req(input$KSEA_data)
    df <- const_df()
    file_name <- input$KSEA_data$name
    new_row <- data.frame(Var = "KSEA_data", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$TreeSVG, {
    req(input$TreeSVG)
    df <- const_df()
    file_name <- input$TreeSVG$name
    new_row <- data.frame(Var = "tree_path", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$upload_meta, {
    req(input$upload_meta)
    df <- const_df()
    file_name <- input$upload_meta$name
    df <- df[df$Var != "Meta1", ]
    new_row <- data.frame(Var = "Meta1", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$upload_meta2, {
    req(input$upload_meta2)
    df <- const_df()
    file_name <- input$upload_meta2$name
    df <- df[df$Var != "Meta2", ]
    new_row <- data.frame(Var = "Meta2", Select = file_name, stringsAsFactors = FALSE)
    const_df(rbind(df, new_row))
  })
  
  observeEvent(input$apply_filter, {
    req(data(), meta(), input$filter_num)
    df <- const_df()
    df <- df[!df$Var %in% c("FilterNum", "FilterCond"), ]
    new_rows <- data.frame(
      Var = c("FilterNum", "FilterCond"),
      Select = c(input$filter_num, input$filterop1),
      stringsAsFactors = FALSE
    )
    const_df(rbind(df, new_rows))
  })
  
  observeEvent(input$apply_filter2, {
    req(data3(), meta2(), input$filter_num2)
    df <- const_df()
    df <- df[!df$Var %in% c("FilterNum2", "FilterCond2"), ]
    new_rows <- data.frame(
      Var = c("FilterNum2", "FilterCond2"),
      Select = c(input$filter_num2, input$filterop2),
      stringsAsFactors = FALSE
    )
    const_df(rbind(df, new_rows))
  })
  
  combined_df <- reactive({
    req(rv$data)
    if (nrow(const_df()) == 0) {
      return(rv$data)
    }
    rbind(rv$data, const_df())
  })
  
  output$log <- downloadHandler(
    filename = function() {
      count <- increment_export_count(counter_file)
      return(paste0("log_", Sys.Date(), "_", count, ".csv"))
    },
    content = function(file) {
      write.csv(combined_df(), file, row.names = FALSE)
    }
  )
  
  output$LogTable <- renderTable({
    df <- combined_df() 
    if (is.data.frame(df)) {  
      xtable(df)  
    } else {
      NULL 
    }
  })
  
  var_map <- list(
    #Color
    Color = "color_palette",
    
    # Meta
    Transformed = "was_transformed",
    Transformed3 = "was_transformed3",
    
    # Imputation
    Imputed = "was_imputed",
    Imputed3 = "was_imputed3",
    
    # Coverage Plot
    CoveragePlotID = "id3",
    CoveragePlotHeader = "header3",
    CoveragePlotLegend = "legend3",
    
    # Missing Value Plot
    MissingValuePlotHeader = "header4",
    MissingValuePlotBin = "missValBin4",
    MissingValueText = "text4",
    MissingValueTextSize = "text_size4",
    
    # Histogram Intensity
    HistogramIntHeader = "header5",
    HistogramIntLegend = "legend5",
    
    # Boxplot Intensity
    BoxplotIntID = "id6",
    BoxplotIntOut = "outliersBX",
    BoxplotIntMean = "meanBX",
    BoxplotIntHeader = "header6",
    BoxplotIntLegend = "legend6",
    
    # COV Plot
    CovPlotOut = "outliersCOV",
    CovPlotHeader = "header7",
    CovPlotLegend = "legend7",
    
    # PCA Plot
    PCAHeader = "header8",
    PCALegend = "legend8",
    
    # CorrPlot
    CorrPlotDisplay = "MeCorr",
    CorrPlotID = "id12",
    
    # Heatmap
    HeatmapID = "id13",
    
    # RT Plot
    RTPlotSytle = "style14",
    RTline = "line14",
    HexbinsRT = "bins14",
    RTHeader = "header14",
    
    # Modifications Plot
    ModPlotID = "id15",
    ModPlotHeader = "header15",
    
    # Missed Cleavage Plot
    MissedCleavID = "id16",
    MissedCleavText = "text16",
    MissedCleavTextSize = "text_size16",
    MissedCleavHeader = "header16"
  )
  
  observeEvent(input$upload_log, {
    req(input$upload_log)
    
    log_df <- read.csv(input$upload_log$datapath, stringsAsFactors = FALSE)
    
    log_version_row <- log_df$Select[log_df$Var == "Version"]
    log_version <- if (length(log_version_row) > 0 && nzchar(log_version_row)) log_version_row else "V0.0.1"
    
    used_version <- "V0.1.1"
    if (log_version != used_version) {
      msg <- paste0("⚠️ Wrong Proteoics Copilot version used: expected ", log_version,
                    ", but used ", used_version, "! Please use the correct version.")
      showNotification(msg, type = "warning", duration = NULL)
    }
    
    available_file_keys <- c("Data1", "Data2", "Data3", "Meta1", "Meta2")
    available_expected_names <- list()
    
    for (key in available_file_keys) {
      if (key %in% log_df$Var) {
        val <- log_df$Select[log_df$Var == key]
        if (!is.na(val) && nzchar(val)) {
          available_expected_names[[key]] <- val
        }
      }
    }
    
    uploaded_names <- list(
      Data1 = if (!is.null(input$file)) input$file$name else NA,
      Data2 = if (!is.null(input$file2)) input$file2$name else NA,
      Data3 = if (!is.null(input$file3)) input$file3$name else NA,
      Meta1 = if (!is.null(input$upload_meta)) input$upload_meta$name else NA,
      Meta2 = if (!is.null(input$upload_meta2)) input$upload_meta2$name else NA
    )
    
    mismatches <- names(available_expected_names)[
      sapply(names(available_expected_names), function(key) {
        expected <- available_expected_names[[key]]
        uploaded <- uploaded_names[[key]]
        is.na(uploaded) || uploaded != expected
      })
    ]
    
    if (length(mismatches) == 0) {
      showNotification("✅ Log file successfully matched with uploaded files.", type = "message")
    } else {
      msg <- paste("❌ Mismatch or missing upload for:", paste(mismatches, collapse = ", "))
      showNotification(msg, type = "error", duration = NULL)
    }
    
    for (i in seq_len(nrow(log_df))) {
      log_var <- log_df$Var[i]
      val <- log_df$Select[i]
      
      if (!is.na(val) && tolower(val) %in% c("true", "false")) {
        if (log_var %in% c("MissingValuePlotBin", "HexbinsRT")) {
          val <- ifelse(tolower(val) == "false", 0, 1)
        } else {
          val <- tolower(val) == "true"
        }
      } else if (!is.na(suppressWarnings(as.numeric(val)))) {
        val <- as.numeric(val)
      } else {
        val <- as.character(val)
      }
      
      reactive_var_name <- var_map[[log_var]]
      
      if (!is.null(reactive_var_name)) {
        if (reactive_var_name %in% c("style14", "color_palette")) {
          updateSelectInput(session, reactive_var_name, selected = val)
        } else if (reactive_var_name %in% c("bins14", "missValBin4", "text_size4", "text_size16")) {
          updateNumericInput(session, reactive_var_name, value = val)
        } else {
          reactive_val <- get(reactive_var_name, envir = environment())
          reactive_val(val)
        }
      }
    }
    
    additional_rows <- log_df[grepl("text\\d+(up|down)", log_df$Var) | grepl("VolcPlot", log_df$Var), ]
    if (nrow(additional_rows) > 0) {
      current_df <- const_df()
      current_df <- current_df[!current_df$Var %in% additional_rows$Var, ]
      updated_df <- rbind(current_df, additional_rows)
      const_df(updated_df)
    }
  })
  
  observeEvent(input$upload_log, {
    req(data(), meta())
    
    log_df <- read.csv(input$upload_log$datapath, stringsAsFactors = FALSE)
    
    transformed_flag1 <- FALSE
    if ("Transformed" %in% log_df$Var) {
      transformed_flag1 <- tolower(log_df$Select[log_df$Var == "Transformed"]) == "false"
    }
    
    if (transformed_flag1) {
      transformed_data(data())
      not_transformed_data(inverseof_log2_transform_data(data(), meta()))
      was_transformed(FALSE)
    } else {
      transformed_data(log2_transform_data(data(), meta()))
      not_transformed_data(data())
      was_transformed(TRUE)
    }
  })
  
  observeEvent(input$upload_log, {
    req(data(), transformed_data(), meta())
    
    log_df <- read.csv(input$upload_log$datapath, stringsAsFactors = FALSE)
    
    current_df <- const_df()
    current_df <- current_df[!current_df$Var %in% c("FilterNum", "FilterCond"), ]
    print(current_df)
    
    if (all(c("FilterNum", "FilterCond") %in% log_df$Var)) {
      print("IWTKMS")
      filter_num_val <- suppressWarnings(as.numeric(log_df$Select[log_df$Var == "FilterNum"]))
      filter_cond_val <- as.character(log_df$Select[log_df$Var == "FilterCond"])
      
      if (!is.na(filter_num_val) && nzchar(filter_cond_val)) {
        new_rows <- data.frame(
          Var = c("FilterNum", "FilterCond"),
          Select = c(filter_num_val, filter_cond_val),
          stringsAsFactors = FALSE
        )
        updated_df <- rbind(current_df, new_rows)
        const_df(updated_df)
        
        filtered_data(filter_data(data(), meta(), num = filter_num_val, filterops = filter_cond_val))
        log2_filtered_data(filter_data(transformed_data(), meta(), num = filter_num_val, filterops = filter_cond_val))
        
        updateNumericInput(session, "filter_num", value = filter_num_val)
        updateSelectInput(session, "filterop1", selected = filter_cond_val)
      } else {
        const_df(current_df)
      }
    } else {
      const_df(current_df)
    }
  })
  
  observeEvent(input$upload_log, {
    
    volcano_data <- reactive({
      if (!is.null(imputed_data())) {
        req(imputed_data())
        return(imputed_data())
      } else {
        if (!is.null(log2_filtered_data())) {
          return(log2_filtered_data())
        } else {
          req(transformed_data())
          return(transformed_data())
        }
      }
    })
    
    req(input$upload_log, pre_volcano_data(), meta())
    
    log_df <- read.csv(input$upload_log$datapath, stringsAsFactors = FALSE)
    
    current_df <- const_df()
    current_df <- current_df[!current_df$Var %in% c("Q-shift", "Adj-std", "Seed"), ]
    print(current_df)
    
    if (all(c("Q-shift", "Adj-std", "Seed") %in% log_df$Var)) {
      q_val <- suppressWarnings(as.numeric(log_df$Select[log_df$Var == "Q-shift"]))
      adj_std_val <- suppressWarnings(as.numeric(log_df$Select[log_df$Var == "Adj-std"]))
      seed_val <- suppressWarnings(as.integer(log_df$Select[log_df$Var == "Seed"]))
      
      if (!is.na(q_val) && !is.na(adj_std_val) && !is.na(seed_val)) {
        new_rows <- data.frame(
          Var = c("Q-shift", "Adj-std", "Seed"),
          Select = c(q_val, adj_std_val, seed_val),
          stringsAsFactors = FALSE
        )
        
        updated_df <- rbind(current_df, new_rows)
        const_df(updated_df)
        
        imputed_data(impute_values(
          pre_volcano_data(), meta(),
          ret = 0,
          q = q_val,
          adj_std = adj_std_val,
          seed = seed_val
        ))
        volcano_data <- reactive({ imputed_data() })
        was_imputed(TRUE)
        
        updateNumericInput(session, "qn1", value = q_val)
        updateNumericInput(session, "adj_stdn1", value = adj_std_val)
        updateNumericInput(session, "seedn1", value = seed_val)
        
      } else {
        const_df(current_df)
      }
    } else {
      const_df(current_df)
    }
  })
  
  observeEvent(input$upload_log, {
    req(data3(), meta2())
    
    log_df <- read.csv(input$upload_log$datapath, stringsAsFactors = FALSE)
    
    transformed_flag3 <- FALSE
    if ("Transformed3" %in% log_df$Var) {
      transformed_flag3 <- tolower(log_df$Select[log_df$Var == "Transformed3"]) == "false"
    }
    
    if (transformed_flag3) {
      transformed_data3(data3())
      not_transformed_data3(inverseof_log2_transform_data(data3(), meta2()))
      was_transformed3(FALSE)
    } else {
      transformed_data3(log2_transform_data(data3(), meta2()))
      not_transformed_data3(data3())
      was_transformed3(TRUE)
    }
  })
  
  observeEvent(input$upload_log, {
    req(data3(), transformed_data3(), meta2())
    
    log_df <- read.csv(input$upload_log$datapath, stringsAsFactors = FALSE)
    
    current_df <- const_df()
    current_df <- current_df[!current_df$Var %in% c("FilterNum2", "FilterCond2"), ]
    
    if (all(c("FilterNum2", "FilterCond2") %in% log_df$Var)) {
      filter_num_val2 <- suppressWarnings(as.numeric(log_df$Select[log_df$Var == "FilterNum2"]))
      filter_cond_val2 <- as.character(log_df$Select[log_df$Var == "FilterCond2"])
      
      new_rows <- data.frame(
        Var = c("FilterNum2", "FilterCond2"),
        Select = c(filter_num_val2, filter_cond_val2),
        stringsAsFactors = FALSE
      )
      updated_df <- rbind(current_df, new_rows)
      const_df(updated_df)
      
      if (!is.na(filter_num_val2) && nzchar(filter_cond_val2)) {
        filtered_data3(filter_data(data3(), meta2(), num = filter_num_val2, filterops = filter_cond_val2))
        log2_filtered_data3(filter_data(transformed_data3(), meta2(), num = filter_num_val2, filterops = filter_cond_val2))
        
        updateNumericInput(session, "filter_num3", value = filter_num_val2)
        updateSelectInput(session, "filterop3", selected = filter_cond_val2)
      }
    } else {
      const_df(current_df)
    }
  })
  
  observeEvent(input$upload_log, {
    req(input$upload_log, volcano_data3(), meta2())
    
    log_df <- read.csv(input$upload_log$datapath, stringsAsFactors = FALSE)
    
    current_df <- const_df()
    current_df <- current_df[!current_df$Var %in% c("Qn2", "AdjStd2", "Seed2"), ]
    
    if (all(c("Qn2", "AdjStd2", "Seed2") %in% log_df$Var)) {
      q_val2 <- suppressWarnings(as.numeric(log_df$Select[log_df$Var == "Qn2"]))
      adj_std_val2 <- suppressWarnings(as.numeric(log_df$Select[log_df$Var == "AdjStd2"]))
      seed_val2 <- suppressWarnings(as.integer(log_df$Select[log_df$Var == "Seed2"]))
      
      if (!is.na(q_val2) && !is.na(adj_std_val2) && !is.na(seed_val2)) {
        new_rows <- data.frame(
          Var = c("Qn2", "AdjStd2", "Seed2"),
          Select = c(q_val2, adj_std_val2, seed_val2),
          stringsAsFactors = FALSE
        )
        
        updated_df <- rbind(current_df, new_rows)
        const_df(updated_df)
        
        imputed_data3(impute_values(
          volcano_data3(), meta2(),
          ret = 0,
          q = q_val2,
          adj_std = adj_std_val2,
          seed = seed_val2
        ))
        volcano_data3 <- reactive({ imputed_data3() })
        was_imputed3(TRUE)
        
        updateNumericInput(session, "qn2", value = q_val2)
        updateNumericInput(session, "adj_stdn2", value = adj_std_val2)
        updateNumericInput(session, "seedn2", value = seed_val2)
        
      } else {
        const_df(current_df)
      }
    } else {
      const_df(current_df)
    }
  })
  
  observeEvent(input$upload_log, {
    
    volcano_data3 <- reactive({
      if (!is.null(imputed_data3())) {
        req(imputed_data3())
        return(imputed_data3())
      } else {
        if (!is.null(log2_filtered_data3())) {
          return(log2_filtered_data3())
        } else {
          req(transformed_data3())
          return(transformed_data3())
        }
      }
    })
    
    req(input$upload_log, pre_volcano_data3(), meta2())
    
    log_df <- read.csv(input$upload_log$datapath, stringsAsFactors = FALSE)
    
    current_df <- const_df()
    current_df <- current_df[!current_df$Var %in% c("Q-shift2", "Adj-std2", "Seed2"), ]
    
    if (all(c("Q-shift2", "Adj-std2", "Seed2") %in% log_df$Var)) {
      q_val2 <- suppressWarnings(as.numeric(log_df$Select[log_df$Var == "Q-shift2"]))
      adj_std_val2 <- suppressWarnings(as.numeric(log_df$Select[log_df$Var == "Adj-std2"]))
      seed_val2 <- suppressWarnings(as.integer(log_df$Select[log_df$Var == "Seed2"]))
      
      if (!is.na(q_val2) && !is.na(adj_std_val2) && !is.na(seed_val2)) {
        new_rows <- data.frame(
          Var = c("Q-shift2", "Adj-std2", "Seed2"),
          Select = c(q_val2, adj_std_val2, seed_val2),
          stringsAsFactors = FALSE
        )
        
        updated_df <- rbind(current_df, new_rows)
        const_df(updated_df)
        
        imputed_data3(impute_values(
          pre_volcano_data3(), meta2(),
          ret = 0,
          q = q_val2,
          adj_std = adj_std_val2,
          seed = seed_val2
        ))
        volcano_data3 <- reactive({ imputed_data3() })
        was_imputed3(TRUE)
        
        updateNumericInput(session, "qn1", value = q_val2)
        updateNumericInput(session, "adj_stdn1", value = adj_std_val2)
        updateNumericInput(session, "seedn1", value = seed_val2)
        
      } else {
        const_df(current_df)
      }
    } else {
      const_df(current_df)
    }
  })
  
  
}

#### ShinyApp logic ####
shinyApp(ui = ui, server = server)