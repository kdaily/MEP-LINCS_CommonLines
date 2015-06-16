library(dplyr)
library(tidyr)
library(data.table)
library(synapseClient)
synapseLogin()

# CCLE Gene Expression
ccleFullGELineNames <- c("MCF7_BREAST", "PC3_PROSTATE", "HPAC_PANCREAS", "A375_SKIN", "A549_LUNG")
ccleFullGEObj <- synGet("syn2318265")
ccleFullGE <- readRDS(getFileLocation(ccleFullGEObj))
ccleFullGELINCS <- ccleFullGE[, ccleFullGELineNames]

# CCLE Gene Expression (maybe processed?)
ccleGELineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")
ccleGEObj <- synGet("syn2293298")
ccleGE <- fread(getFileLocation(ccleGEObj), data.table=FALSE)
ccleGELINCS <- ccleGE[, ccleGELineNames]

# CCLE Copy Number
ccleCNALineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")
ccleCNAObj <- synGet("syn2293299")
ccleCNA <- fread(getFileLocation(ccleCNAObj), data.table=FALSE)
ccleCNALINCS <- ccleCNA[, ccleCNALineNames]

# CCLE Drug EC50
# Data is cell lines (rows) x drugs (columns)
ccleEC50LineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")
ccleEC50Obj <- synGet("syn2293304")
ccleEC50 <- fread(getFileLocation(ccleEC50Obj), data.table=FALSE)
rownames(ccleEC50) <- ccleEC50$V1
ccleEC50$V1 <- NULL
ccleEC50LINCS <- t(ccleEC50[ccleEC50LineNames, ])

# CCLE Drug IC50
# Data is cell lines (rows) x drugs (columns)
ccleIC50LineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")
ccleIC50Obj <- synGet("syn2293303")
ccleIC50 <- fread(getFileLocation(ccleIC50Obj), data.table=FALSE)
rownames(ccleIC50) <- ccleIC50$V1
ccleIC50$V1 <- NULL
ccleIC50LINCS <- t(ccleIC50[ccleIC50LineNames, ])


# Broad CTD2
ctd2LineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")

# AUC
ctd2AUCObj <- synGet("syn4260138")
ctd2AUC <- fread(getFileLocation(ctd2AUCObj), data.table=FALSE)
ctd2AUCLINCS <- ctd2AUC %>% filter(cell_line_name %in% ctd2LineNames)

# Raw viability
ctd2RawViabilityObj <- synGet("syn4260136")
ctd2RawViability <- fread(getFileLocation(ctd2RawViabilityObj), data.table=FALSE)
ctd2RawViabilityLINCS <- ctd2RawViability %>% filter(cell_line_name %in% ctd2LineNames)
