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

get_current_script_path <- function() {
  tryCatch({
    normalizePath(attr(attr(parent.frame(), "srcref"), "srcfile")$filename)
  }, error = function(e1) {
    tryCatch({
      normalizePath(rstudioapi::getSourceEditorContext()$path)
    }, error = function(e2) {
      tryCatch({
        args <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("^--file=", args, value = TRUE)
        if (length(file_arg) > 0) {
          normalizePath(sub("^--file=", "", file_arg))
        } else {
          stop("Cannot determine the path of the running script")
        }
      }, error = function(e3) {
        stop("Cannot determine the path of the running script")
      })
    })
  })
}

cat("Current script path:\n", get_current_script_path(), "\n")