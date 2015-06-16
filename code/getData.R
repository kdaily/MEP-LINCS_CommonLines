library(dplyr)
library(tidyr)
library(synapseClient)
synapseLogin()

ccleFullGEObj <- synGet("syn2318265")
ccleFullGE <- readRDS(getFileLocation(ccleFullGEObj))

ccleFullGELineNames <- c("MCF7_BREAST", "PC3_PROSTATE", "HPAC_PANCREAS", "A375_SKIN", "A549_LUNG")
ccleLINCSGE <- ccleFullGE[, ccleFullGELineNames]
