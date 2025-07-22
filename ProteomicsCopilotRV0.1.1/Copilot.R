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

copilot_path <- tryCatch(
  normalizePath(sys.frame(1)$ofile),
  error = function(e) normalizePath(rstudioapi::getSourceEditorContext()$path)
)

app_dir <- dirname(copilot_path)

cat("Copilot path:\n", copilot_path, "\n")
cat("App directory contents:\n")
print(list.files(app_dir))

if (!file.exists(file.path(app_dir, "app.R"))) {
  stop("No app.R found in the same folder as Copilot.R")
}

setwd(app_dir)
shiny::runApp(app_dir)

#setwd("C:/Users/jonas/OneDrive/Desktop/Vis_phos/Copilot Shiny2")
#shinylive::export(appdir = "web app", destdir = "docs")
#httpuv::runStaticServer("docs/", port=8008)
