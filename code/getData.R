library(dplyr)
library(tidyr)
library(data.table)
library(reshape2)
library(synapseClient)
library(rGithubClient)
synapseLogin()

repo <- getRepo("kdaily/MEP-LINCS_Common_Lines", ref="branch", refName="master")
script <- getPermlink(repo, "code/getData.R")

projectId <- "syn4259323"
proj <- synGet(projectId)

# CCLE
ccleLineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")

# CCLE Gene Expression (maybe processed?)
id <- "syn2293298"
ccleGEObj <- synGet(id)
ccleGE <- fread(getFileLocation(ccleGEObj), data.table=FALSE)
ccleGELINCS <- ccleGE[, ccleLineNames]
write.csv(ccleGELINCS, file="CCLE_GE_LINCS.csv")
f <- File("CCLE_GE_LINCS.csv", parentId=proj)
act <- Activity(used=id, executed=thisScript)

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

# JW Gray
grayLineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")
grayExonObj <- synGet("syn2346647")
grayExon <- fread(getFileLocation(grayExonObj), data.table=FALSE)
grayExonLINCS <- grayExon[, c("GeneSymbol", intersect(grayLineNames, colnames(grayExon)))]

grayRPPAObj <- synGet("syn2347012")
grayRPPA <- fread(getFileLocation(grayRPPAObj), data.table=FALSE)
grayRPPAMeta <- grayRPPA[1, ]
t(as.data.frame(grayRPPA[-1, ]))
grayRPPALINCS <- grayRPPA[, c("GeneSymbol", intersect(grayLineNames, colnames(grayRPPA)))]

graySNP6Obj <- synGet("syn2347009")
graySNP6 <- fread(getFileLocation(graySNP6Obj), data.table=FALSE)
graySNP6LINCS <- graySNP6[, c("chrom", "start", "end", intersect(grayLineNames, colnames(graySNP6)))]

grayWesternObj <- synGet("syn2347011")
grayWestern <- fread(getFileLocation(grayWesternObj), data.table=FALSE)
grayWesternLINCS <- grayWestern[, c("Protein", intersect(grayLineNames, colnames(grayWestern)))]

grayDblObj <- synGet("syn2347014")
grayDbl <- fread(getFileLocation(grayDblObj), data.table=FALSE)
grayDblLINCS <- grayDbl %>% filter(CellLineName %in% grayLineNames)

grayRNASeqMatObj <- synGet("syn2347004")
grayRNASeqMat <- fread(getFileLocation(grayRNASeqMatObj), data.table=FALSE)

grayRNASeqMatLINCS <- grayRNASeqMat %>% select(Gene_ID, FID, Seq_Name, EnsEMBL_Gene_ID, 
                                               one_of(intersect(grayLineNames, colnames(grayWestern))))
