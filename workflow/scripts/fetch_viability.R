library(PharmacoGx)
library(tidyr)
library(dplyr)
library(tibble)
library(yaml)

project.root <- Sys.getenv("PIXI_PROJECT_ROOT")
config <- read_yaml(file.path(Sys.getenv("CONFIG"), "downloads.yaml"))

config <- config$VIABILITY

write.dir <- file.path(project.root, config$dir)


raw.dir <- file.path(project.root, config$rawdir, "PSETS")
pset.name <- config$screen

# for(ps.name in c("NCI60","CTRPv2")){
#     print(paste0("Working on ",ps.name))
dir.create(raw.dir, showWarnings = FALSE, recursive = TRUE)
dir.create(write.dir, showWarnings = FALSE, recursive = TRUE)
ps <- downloadPSet(pset.name, saveDir = raw.dir)
ps <- updateObject(ps)


treatment.info <- treatmentInfo(ps) |> select("treatmentid", "PubChem.CID")
aucs <- summarizeSensitivityProfiles(ps, sensitivity.measure = 'aac_recomputed')
df <- data.frame(aucs)
df <- rownames_to_column(df, var = "treatmentid")
df <- merge(treatment.info, df, by = "treatmentid")
write.csv(df, file.path(write.dir, paste0(pset.name, ".csv")))
