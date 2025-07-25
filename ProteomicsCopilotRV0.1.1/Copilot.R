#V0.1.1
all_pkgs <- c(
  "shiny", "readr", "readxl", "dplyr", "reshape2", "ggplot2",
  "stringr", "plotly", "tidyverse", "gprofiler2", "corrplot",
  "pheatmap", "reticulate", "ggeasy", "zip", "DT", "rsvg",
  "protr", "r3dmol", "KSEAapp", "xtable", "tools", "png",
  "UniprotR", "protti", "shinycssloaders",
  "Biostrings", "GenomicAlignments", "rstudioapi", "seqinr",
  "openxlsx", "this.path"
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
library(this.path)
setwd(this.path::here())
runApp()