library(dplyr)
library(tidyr)
library(data.table)
library(reshape2)
library(synapseClient)
library(rGithubClient)
synapseLogin()

repo <- getRepo("kdaily/MEP-LINCS_CommonLines", ref="branch", refName="master")
thisScript <- getPermlink(repo, "code/getData.R")

projectId <- "syn4259323"
proj <- synGet(projectId)

# CCLE
ccleLineNames <- c("MCF7", "PC3", "YAPC", "A375", "A549")

# CCLE Gene Expression (maybe processed?)
id <- "syn2293298"
fn <- "CCLE_GE_LINCS.csv"
ccleGEObj <- synGet(id)

# Duplicate columns (though diff cell lines) - get rid of them
ccleGE <- fread(getFileLocation(ccleGEObj), drop=c("TT", "NCIH292"), data.table=FALSE)
ccleGELINCS <- ccleGE %>% 
  select(V1, one_of(ccleLineNames)) %>%
  dplyr::rename(GeneSymbol=V1)

write.csv(ccleGELINCS, file=fn)
f <- File(fn, name="CCLE Gene Expression", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="CCLE", dataType="mRNA", fileType="genomicMatrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

# CCLE Copy Number
id <- "syn2293299"
fn <- "CCLE_CNA_LINCS.csv"
ccleCNAObj <- synGet(id)
ccleCNA <- fread(getFileLocation(ccleCNAObj), data.table=FALSE)
ccleCNALINCS <- ccleCNA[, c("V1", ccleLineNames)] %>%
  rename(GeneSymbol=V1)
write.csv(ccleCNALINCS, file=fn)
f <- File(fn, name="CCLE CNA", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="CCLE", dataType="DNA", 
                             dataSubType="CNA", fileType="genomicMatrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

# CCLE Drug EC50
# Data is cell lines (rows) x drugs (columns)
id <- "syn2293304"
fn <- "CCLE_EC50_LINCS.csv"
ccleEC50Obj <- synGet(id)
ccleEC50 <- fread(getFileLocation(ccleEC50Obj), data.table=FALSE)
ccleEC50LINCS <- ccleEC50 %>%
  rename(cell_line_name=V1) %>% filter(cell_line_name %in% ccleLineNames) %>%
  melt(id.vars='cell_line_name') %>%
  dcast(variable ~ cell_line_name, value.var = 'value') %>%
  rename(drug=variable)
write.csv(ccleEC50LINCS, file=fn)
f <- File(fn, name="CCLE EC50", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="CCLE", dataType="drug",
                             dataSubType="EC50", fileType="matrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

# CCLE Drug IC50
# Data is cell lines (rows) x drugs (columns)
id <- "syn2293303"
fn <- "CCLE_IC50_LINCS.csv"
ccleIC50Obj <- synGet(id)
ccleIC50 <- fread(getFileLocation(ccleIC50Obj), data.table=FALSE)
ccleIC50LINCS <- ccleIC50 %>%
  rename(cell_line_name=V1) %>% filter(cell_line_name %in% ccleLineNames) %>%
  melt(id.vars='cell_line_name') %>%
  dcast(variable ~ cell_line_name, value.var = 'value') %>%
  rename(drug=variable)
write.csv(ccleIC50LINCS, file=fn)
f <- File(fn, name="CCLE IC50", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="CCLE", dataType="drug",
                             dataSubType="IC50", fileType="matrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

# Broad CTD2
ctd2LineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")

# AUC
id <- "syn4260138"
fn <- "CTD2_AUC_LINCS.csv"
ctd2AUCObj <- synGet(id)
ctd2AUC <- fread(getFileLocation(ctd2AUCObj), data.table=FALSE)
ctd2AUCLINCS <- ctd2AUC %>% filter(cell_line_name %in% ctd2LineNames)
write.csv(ctd2AUCLINCS, file=fn, row.names=FALSE)
f <- File(fn, name="CTD2 AUC", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="CTD2", dataType="drug",
                             dataSubType="AUC", fileType="matrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

# Raw viability
id <- "syn4260136"
fn <- "CTD2_RawViability_LINCS.csv"
ctd2RawViabilityObj <- synGet(id)
ctd2RawViability <- fread(getFileLocation(ctd2RawViabilityObj), data.table=FALSE)
ctd2RawViabilityLINCS <- ctd2RawViability %>% filter(cell_line_name %in% ctd2LineNames)
write.csv(ctd2RawViabilityLINCS, file=fn, row.names=FALSE)
f <- File(fn, name="CTD2 Raw Viability", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="CTD2", dataType="drug",
                             dataSubType="viability", fileType="matrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

# Avg pct viability
id <- "syn4260137"
fn <- "CTD2_AvgPctViability_LINCS.csv"
ctd2PctViabilityObj <- synGet(id)
ctd2PctViability <- fread(getFileLocation(ctd2PctViabilityObj), data.table=FALSE)
ctd2PctViabilityLINCS <- ctd2PctViability %>% filter(cell_line_name %in% ctd2LineNames)
write.csv(ctd2PctViabilityLINCS, file=fn, row.names=FALSE)
f <- File(fn, name="CTD2 Average Percent Viability", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="CTD2", dataType="drug",
                             dataSubType="viability", fileType="matrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

# Human Protein Atlas
hpaLineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")

id <- "syn4487737"
hpaObj <- synGet(id)
hpa <- fread(getFileLocation(hpaObj), data.table=FALSE)
hpaLINCS <- hpa %>% filter(Sample %in% hpaLineNames)

fn <- "HPA_FPKM_LINCS.csv"
hpaLINCSFPKM <- hpaLINCS %>% dcast(Gene ~ Sample, value.var="Value")
write.csv(hpaLINCSFPKM, file=fn, row.names=FALSE)
f <- File(fn, name="HPA RNASeq FPKM", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="HPA", dataType="mRNA",
                             dataSubType="FPKM", fileType="genomicMatrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

fn <- "HPA_Abundance_LINCS.csv"
hpaLINCSAbundance <- hpaLINCS %>% dcast(Gene ~ Sample, value.var="Abundance")
write.csv(hpaLINCSAbundance, file=fn, row.names=FALSE)
f <- File(fn, name="HPA RNASeq Abundance", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="HPA", dataType="mRNA",
                             dataSubType="Abundance", fileType="genomicMatrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

# JW Gray
grayLineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")

id <- "syn2346647"
fn <- "JWGray_Exon_LINCS.csv"
grayExonObj <- synGet(id)
grayExon <- fread(getFileLocation(grayExonObj), data.table=FALSE)
grayExonLINCS <- grayExon[, c("GeneSymbol", intersect(grayLineNames, colnames(grayExon)))]
write.csv(grayExonLINCS, file=fn, row.names=FALSE)
f <- File(fn, name="JW Gray Exon", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="JWGray", dataType="mRNA",
                             dataSubType="expr", fileType="genomicMatrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

id <- "syn2347012"
fn <- "JWGray_RPPA_LINCS.csv"
grayRPPAObj <- synGet(id)
grayRPPA <- fread(getFileLocation(grayRPPAObj), data.table=FALSE)

grayRPPALINCS <- grayRPPA[, c("Sample", "FullyValidated", intersect(grayLineNames, colnames(grayRPPA)))]
write.csv(grayRPPALINCS, file=fn, row.names=FALSE)
f <- File(fn, name="JW Gray RPPA", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="JWGray", dataType="protein",
                             dataSubType="RPPA", fileType="genomicMatrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

id <- "syn2347009"
fn <- "JWGray_SNP6_LINCS.csv"
graySNP6Obj <- synGet(id)
graySNP6 <- fread(getFileLocation(graySNP6Obj), data.table=FALSE)
graySNP6LINCS <- graySNP6[, c("chrom", "start", "end", intersect(grayLineNames, colnames(graySNP6)))]
write.csv(graySNP6LINCS, file=fn, row.names=FALSE)
f <- File(fn, name="JW Gray SNP6", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="JWGray", dataType="DNA",
                             dataSubType="CNA", fileType="genomicMatrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

id <- "syn2347011"
fn <- "JWGray_Western_LINCS.csv"
grayWesternObj <- synGet(id)
grayWestern <- fread(getFileLocation(grayWesternObj), data.table=FALSE)
grayWesternLINCS <- grayWestern[, c("Protein", intersect(grayLineNames, colnames(grayWestern)))]
write.csv(grayWesternLINCS, file=fn, row.names=FALSE)
f <- File(fn, name="JW Gray Western Blot", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="JWGray", dataType="protein",
                             dataSubType="Western", fileType="matrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

id <- "syn2347014"
fn <- "JWGray_DblTime_LINCS.csv"
grayDblObj <- synGet(id)
grayDbl <- fread(getFileLocation(grayDblObj), data.table=FALSE)
grayDblLINCS <- grayDbl %>% filter(CellLineName %in% grayLineNames)
write.csv(grayDblLINCS, file=fn, row.names=FALSE)
f <- File(fn, name="JW Gray Doubling Time", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="JWGray", dataType="cell_line",
                             dataSubType="DoublingTime", fileType="matrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

id <- "syn2347004"
fn <- "JWGray_RNASeq_LINCS.csv"
grayRNASeqMatObj <- synGet(id)
grayRNASeqMat <- fread(getFileLocation(grayRNASeqMatObj), data.table=FALSE)

grayRNASeqMatLINCS <- grayRNASeqMat %>% select(Gene_ID, FID, Seq_Name, EnsEMBL_Gene_ID, 
                                               one_of(intersect(grayLineNames, colnames(grayRNASeqMat))))
write.csv(grayRNASeqMatLINCS, file=fn, row.names=FALSE)
f <- File(fn, name="JW Gray RNASeq Expression", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="JWGray", dataType="mRNA",
                             dataSubType="",
                             fileType="genomicMatrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

# Sanger
sangerLineNames <- c("MCF7", "PC3", "HPAC", "A375", "A549")

id <- "syn4513928"
sangerSensObj <- synGet(id)
sangerSense <- fread(getFileLocation(sangerSensObj), data.table=FALSE)

sangerSenseLINCS <- sangerSense %>%
  filter(Cell.Line %in% sangerLineNames)

sangerSenseIC50LINCS <- sangerSenseLINCS %>%
  select(Cell.Line, Cosmic_ID, Cancer.Type, Tissue, ends_with("IC_50"))
colnames(sangerSenseIC50LINCS) <- gsub("_IC_50", "", colnames(sangerSenseIC50LINCS))

fn <- "Sanger_IC50_LINCS.csv"
write.csv(sangerSenseIC50LINCS, file=fn, row.names=FALSE)
f <- File(fn, name="Sanger IC50", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="Sanger", dataType="drug",
                             dataSubType="IC50",
                             fileType="matrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)

sangerSenseAUCLINCS <- sangerSenseLINCS %>%
  select(Cell.Line, Cosmic_ID, Cancer.Type, Tissue, ends_with("AUC"))
colnames(sangerSenseAUCLINCS) <- gsub("_AUC", "", colnames(sangerSenseAUCLINCS))

fn <- "Sanger_AUC_LINCS.csv"
write.csv(sangerSenseAUCLINCS, file=fn, row.names=FALSE)
f <- File(fn, name="Sanger AUC", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="Sanger", dataType="drug",
                             dataSubType="AUC",
                             fileType="matrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)
