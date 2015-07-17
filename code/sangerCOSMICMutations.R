# Getting the file needed from Synapse in this script requires Rsftp
# devtools::install_github("Sage-Bionetworks/Rsftp")

library(dplyr)
library(data.table)
library(synapseClient)
synapseLogin()

projectId <- "syn4259323"
proj <- synGet(projectId)

sangerLineNames <- c("MCF7", "PC3", "YAPC", "A375", "A549")

# Must gunzip this file!
# read.table cannot read it b/c of uneven rows, use fread
# Requires SFTP credentials from Sanger
sangerMutObj <- synGet("syn4633294", downloadLocation='/home/kdaily/data/Sanger/')

system(sprintf("gunzip %s", getFileLocation(sangerMutObj)), 
       intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE, wait=TRUE)
fileLoc <- gsub(".gz", "", getFileLocation(sangerMutObj))

sangerMut <- fread(fileLoc, data.table=FALSE)

sangerMutLINCS <- sangerMut %>% filter(`Sample name` %in% sangerLineNames)
