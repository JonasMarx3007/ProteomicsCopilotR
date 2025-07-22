library(git2r)

repo_url <- "https://github.com/JonasMarx3007/ProteomicsCopilotR.git"
local_path <- "your_path"
temp_path <- tempfile("copilot_")

extract_version_parts <- function(name) {
  ver_str <- sub("^ProteomicsCopilotRV", "", name)
  as.numeric(strsplit(ver_str, "\\.")[[1]])
}

get_version_folder <- function(path) {
  dirs <- list.dirs(path, recursive = FALSE, full.names = FALSE)
  version_dir <- dirs[grepl("^ProteomicsCopilotRV", dirs)]
  if (length(version_dir) == 0) return(NULL)
  return(version_dir[1])
}

read_log_counter <- function(path) {
  log_path <- file.path(path, "www", "ressources", "logcounter.txt")
  if (file.exists(log_path)) {
    tryCatch(as.integer(readLines(log_path, warn = FALSE)[1]), error = function(e) NULL)
  } else NULL
}

write_log_counter <- function(path, value) {
  log_path <- file.path(path, "www", "ressources", "logcounter.txt")
  dir.create(dirname(log_path), recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(value), log_path)
}

local_version_folder <- NULL
local_counter <- NULL

if (dir.exists(local_path)) {
  local_version_folder <- get_version_folder(local_path)
  if (!is.null(local_version_folder)) {
    message("Local version detected: ", local_version_folder)
    local_counter <- read_log_counter(file.path(local_path, local_version_folder))
    if (!is.null(local_counter)) {
      message("Extracted local logcounter value: ", local_counter)
    }
  }
}

message("Fetching newest remote version into temporary folder...")
message("Temporary path used: ", temp_path)
clone(repo_url, temp_path, progress = FALSE)
remote_version_folder <- get_version_folder(temp_path)
message("Remote version folder detected: ", remote_version_folder)

safe_delete <- function(path) {
  try({
    for (i in 1:5) {
      unlink(path, recursive = TRUE, force = TRUE)
      Sys.sleep(1)
      if (!dir.exists(path)) break
    }
  }, silent = TRUE)
  if (dir.exists(path)) {
    stop("Could not delete locked folder: ", path)
  }
}

install_new_version <- function() {
  if (startsWith(normalizePath(getwd(), winslash = "/"), normalizePath(local_path, winslash = "/"))) {
    setwd(dirname(local_path))
  }
  if (dir.exists(local_path)) {
    safe_delete(local_path)
  }
  message("Installing new version...")
  clone(repo_url, local_path)
  if (!is.null(local_counter)) {
    write_log_counter(file.path(local_path, remote_version_folder), local_counter)
    message("Reinserted logcounter value into new installation: ", local_counter)
  }
}

if (is.null(local_version_folder)) {
  install_new_version()
} else {
  local_ver <- extract_version_parts(local_version_folder)
  remote_ver <- extract_version_parts(remote_version_folder)
  if (!identical(local_ver, remote_ver)) {
    if (!all(local_ver[1:2] == remote_ver[1:2])) {
      prompt <- readline("A major or minor update is available. Install? (y/n): ")
      if (tolower(prompt) == "y") {
        install_new_version()
      }
    } else {
      message("Patch version difference detected. Pulling latest changes...")
      repo <- repository(local_path)
      fetch(repo, name = "origin")
      current_branch <- repository_head(repo)$name
      merge(repo, paste0("origin/", current_branch))
      message("Pulled latest changes.")
    }
  } else {
    message("Versions match. Pulling to ensure latest sync...")
    repo <- repository(local_path)
    fetch(repo, name = "origin")
    current_branch <- repository_head(repo)$name
    merge(repo, paste0("origin/", current_branch))
    message("Pulled latest changes.")
  }
}

setwd(local_path)

all_files <- list.files(path = ".", pattern = "Copilot\\.R$", recursive = TRUE)

if (length(all_files) == 0) {
  stop("No 'Copilot.R' found!")
} else if (length(all_files) > 1) {
  warning("More than one 'Copilot.R' found!")
}

copilot_path <- all_files[1]
full_path <- normalizePath(copilot_path, winslash = "/")
message("Sourcing Copilot.R von: ", full_path)
safe_delete(temp_path)
source(copilot_path)
