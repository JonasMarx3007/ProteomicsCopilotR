all_pkgs <- c(
  "shiny", "readr", "readxl", "dplyr", "reshape2", "ggplot2",
  "stringr", "plotly", "tidyverse", "gprofiler2", "corrplot",
  "pheatmap", "reticulate", "ggeasy", "zip", "DT", "rsvg",
  "protr", "r3dmol", "KSEAapp", "xtable", "tools", "png",
  "UniprotR", "protti", "shinycssloaders",
  "Biostrings", "GenomicAlignments", "rstudioapi", "seqinr",
  "openxlsx"
)

bioc_pkgs <- c("Biostrings", "GenomicAlignments", "protti")
cran_pkgs <- setdiff(all_pkgs, bioc_pkgs)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
  library(pkg, character.only = TRUE)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
  library(pkg, character.only = TRUE)
}

library(shiny)

get_copilot_path <- function() {
  tryCatch({
    normalizePath(attr(attr(parent.frame(), "srcref"), "srcfile")$filename)
  }, error = function(e1) {
    tryCatch({
      normalizePath(rstudioapi::getSourceEditorContext()$path)
    }, error = function(e2) {
      tryCatch({
        script_args <- commandArgs(trailingOnly = FALSE)
        script_file <- sub("--file=", "", script_args[grep("--file=", script_args)])
        normalizePath(script_file)
      }, error = function(e3) {
        stop("Could not determine path to Copilot.R")
      })
    })
  })
}

copilot_path <- get_copilot_path()
app_dir <- dirname(copilot_path)
cat("Copilot path:\n", copilot_path, "\n")
cat("App directory contents:\n")
print(list.files(app_dir))

if (!file.exists(file.path(app_dir, "app.R"))) {
  stop("No app.R found in the same folder as Copilot.R")
}

setwd(app_dir)
shiny::runApp(app_dir)
