library(dplyr)
library(tidyr)
library(data.table)
library(reshape2)
library(synapseClient)
library(rGithubClient)
synapseLogin()

repo <- getRepo("kdaily/MEP-LINCS_CommonLines", ref="branch", refName="sangermut")
thisScript <- getPermlink(repo, "code/sangerCOSMICMutations.R")

projectId <- "syn4259323"
proj <- synGet(projectId)

# Sanger
sangerLineNames <- c("MCF7", "PC3", "YAPC", "A375", "A549")

# Must gunzip this file!
# read.table cannot read it b/c of uneven rows, use fread
sangerMutObj <- synGet("syn4633294", downloadLocation='/home/kdaily/data/Sanger/')
system(sprintf("gunzip %s", getFileLocation(sangerMutObj)), 
       intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE, wait=TRUE)

fileLoc <- gsub(".gz", "", getFileLocation(sangerMutObj))
sangerMut <- fread(fileLoc, data.table=FALSE)

sangerMutLINCS <- sangerMut %>% filter(`Sample name` %in% sangerLineNames)
