library(git2r)

repo_url <- "https://github.com/JonasMarx3007/ProteomicsCopilotR.git"
local_path <- "C:/Users/Jonas Marx/Desktop/CopilotNewestVersion"

if (!dir.exists(local_path) || !file.exists(file.path(local_path, ".git", "config"))) {
  message("Cloning repository...")
  clone(repo_url, local_path)
} else {
  message("Repository already exists. Pulling latest changes...")
  repo <- repository(local_path)
  pull(repo)
}