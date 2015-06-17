library(dplyr)
library(tidyr)
library(data.table)
library(reshape2)
library(synapseClient)
synapseLogin()

# CCLE
ccleLineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")

# CCLE Gene Expression (maybe processed?)
ccleGEObj <- synGet("syn2293298")
ccleGE <- fread(getFileLocation(ccleGEObj), data.table=FALSE)
ccleGELINCS <- ccleGE[, ccleLineNames]

# CCLE Copy Number
ccleCNAObj <- synGet("syn2293299")
ccleCNA <- fread(getFileLocation(ccleCNAObj), data.table=FALSE)
ccleCNALINCS <- ccleCNA[, ccleLineNames]

# CCLE Drug EC50
# Data is cell lines (rows) x drugs (columns)
ccleEC50Obj <- synGet("syn2293304")
ccleEC50LINCS <- fread(getFileLocation(ccleEC50Obj), data.table=FALSE) %>%
  rename(cell_line_name=V1) %>% filter(cell_line_name %in% ccleLineNames) %>%
  melt(id.vars='cell_line_name') %>%
  dcast(variable ~ cell_line_name, value.var = 'value') %>%
  rename(drug=variable)

# CCLE Drug IC50
# Data is cell lines (rows) x drugs (columns)
ccleIC50Obj <- synGet("syn2293303")
ccleIC50LINCS <- fread(getFileLocation(ccleIC50Obj), data.table=FALSE) %>%
  rename(cell_line_name=V1) %>% filter(cell_line_name %in% ccleLineNames) %>%
  melt(id.vars='cell_line_name') %>%
  dcast(variable ~ cell_line_name, value.var = 'value') %>%
  rename(drug=variable)

# Broad CTD2
ctd2LineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")

# AUC
ctd2AUCObj <- synGet("syn4260138")
ctd2AUC <- fread(getFileLocation(ctd2AUCObj), data.table=FALSE)
ctd2AUCLINCS <- ctd2AUC %>% filter(cell_line_name %in% ctd2LineNames)
write.csv(ctd2AUCLINCS, file="CTD2_AUC_LINCS.csv", row.names=FALSE)

# Raw viability
ctd2RawViabilityObj <- synGet("syn4260136")
ctd2RawViability <- fread(getFileLocation(ctd2RawViabilityObj), data.table=FALSE)
ctd2RawViabilityLINCS <- ctd2RawViability %>% filter(cell_line_name %in% ctd2LineNames)

# Avg pct viability
ctd2PctViabilityObj <- synGet("syn4260137")
ctd2PctViability <- fread(getFileLocation(ctd2PctViabilityObj), data.table=FALSE)
ctd2PctViabilityLINCS <- ctd2PctViability %>% filter(cell_line_name %in% ctd2LineNames)

# Human Protein Atlas
hpaLineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")
hpaObj <- synGet("syn4487737")
hpa <- fread(getFileLocation(hpaObj), data.table=FALSE)
hpaLINCS <- hpa %>% filter(Sample %in% hpaLineNames)

hpaLINCSMatrix <- hpaLINCS %>% dcast(Gene ~ Sample, value.var="Value")
hpaLINCSAbundance <- hpaLINCS %>% dcast(Gene ~ Sample, value.var="Abundance")
