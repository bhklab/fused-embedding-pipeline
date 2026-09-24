options(repos = c(CRAN = "https://cloud.r-project.org"))
library.path <- Sys.getenv("R_LIBS_USER")
dir.create(library.path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library.path, .libPaths()))

remotes::install_github(
  "bhklab/PharmacoGx@devel",
  lib = library.path,
  dependencies = c("Depends", "Imports", "LinkingTo"),
  upgrade = "never",
  build_vignettes = FALSE
)
