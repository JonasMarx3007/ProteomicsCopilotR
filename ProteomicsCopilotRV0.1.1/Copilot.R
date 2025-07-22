#V0.1.1
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
this_file <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    normalizePath(sub(needle, "", cmdArgs[match]))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    normalizePath(sys.frames()[[1]]$ofile)
  } else {
    tryCatch(
      normalizePath(rstudioapi::getActiveDocumentContext()$path),
      error = function(e) stop("Cannot determine script path")
    )
  }
}

setwd(dirname(this_file()))
runApp()

#setwd("C:/Users/jonas/OneDrive/Desktop/Vis_phos/Copilot Shiny2")
#shinylive::export(appdir = "web app", destdir = "docs")
#httpuv::runStaticServer("docs/", port=8008)
