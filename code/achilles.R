library(dplyr)
library(tidyr)
library(data.table)
library(reshape2)
library(synapseClient)
library(rGithubClient)
synapseLogin()

# repo <- getRepo("kdaily/MEP-LINCS_CommonLines", ref="branch", refName="master")
# thisScript <- getPermlink(repo, "code/achilles.R")

projectId <- "syn4259323"
proj <- synGet(projectId)

availLines <- c() 

achillesLineNames <- c("MCF7_BREAST", "A549_LUNG")

id <- "syn2293298"
fn <- "CCLE_GE_LINCS.csv"
achillesGEObj <- synGet(id)

