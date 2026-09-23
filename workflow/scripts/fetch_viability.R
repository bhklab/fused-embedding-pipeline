library(PharmacoGx)
library(tidyr)
library(dplyr)
library(tibble)
library(yaml)

config<- read_yaml("../../config/config.yaml")

config <- config$VIABILITY

write.dir <- config$dir


raw.dir <- paste0(config$rawdir,"PSETS")
pset.name <- config$screen

# for(ps.name in c("NCI60","CTRPv2")){
#     print(paste0("Working on ",ps.name))
dir.create(raw.dir,showWarnings = FALSE)
ps <- downloadPSet(pset.name,saveDir = raw.dir)
ps <- updateObject(ps)



treatment.info <- treatmentInfo(ps)%>% select("treatmentid", "PubChem.CID")
aucs <- summarizeSensitivityProfiles(ps,sensitivity.measure='aac_recomputed')
df <- data.frame(aucs)
df <- rownames_to_column(df,var="treatmentid")
df <- merge(treatment.info,df, by="treatmentid")
write.csv(df,paste0(write.dir,pset.name,".csv"))
