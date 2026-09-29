options(repos = BiocManager::repositories())
library.path <- Sys.getenv("R_LIBS_USER")
dir.create(library.path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library.path, .libPaths()))

remotes::install_github(
  paste0("bhklab/PharmacoGx@", Sys.getenv("PHARMACOGX_REF", "devel")),
  lib = library.path,
  dependencies = c("Depends", "Imports", "LinkingTo"),
  upgrade = "never",
  build_vignettes = FALSE
)
