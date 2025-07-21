library(git2r)

repo_url <- "https://github.com/JonasMarx3007/ProteomicsCopilotR.git"
local_path <- "C:/Users/Jonas Marx/Desktop/Productive/NewCopilotInstall"

if (!dir.exists(local_path) || !file.exists(file.path(local_path, ".git", "config"))) {
  message("Cloning repository...")
  clone(repo_url, local_path)
} else {
  message("Repository already exists. Pulling latest changes...")
  repo <- repository(local_path)
  
  fetch(repo, name = "origin")
  
  current_branch <- repository_head(repo)$name
  merge(repo, paste0("origin/", current_branch))
}

setwd(local_path)

all_files <- list.files(path = ".", pattern = "Copilot\\.R$", recursive = TRUE)

if (length(all_files) == 0) {
  stop("Datei 'Copilot.R' wurde im Repository nicht gefunden.")
} else if (length(all_files) > 1) {
  warning("Mehrere Dateien namens 'Copilot.R' gefunden. Nimm die erste.")
}

copilot_path <- all_files[1]
full_path <- normalizePath(copilot_path, winslash = "/")
message("Sourcing Copilot.R von: ", full_path)
source(copilot_path)

