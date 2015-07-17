library(dplyr)
library(tidyr)
library(data.table)
library(reshape2)
library(synapseClient)
library(rGithubClient)
synapseLogin()

repo <- getRepo("kdaily/MEP-LINCS_CommonLines", ref="branch", refName="achilles")
thisScript <- getPermlink(repo, "code/achilles.R")

projectId <- "syn4259323"
proj <- synGet(projectId)

achillesLineNames <- c("MCF7_BREAST", "A549_LUNG")

id <- "syn4598414"
fn <- "Achilles_Gs_LINCS.csv"
achillesGsObj <- synGet(id)

achillesGs <- fread(getFileLocation(achillesGsObj), skip=2, data.table=FALSE,
                    select=c("Name", "Description", achillesLineNames)) %>%
  rename(MCF7=MCF7_BREAST, A549=A549_LUNG)

write.csv(achillesGs, file=fn)
f <- File(fn, name="Achilles QC rnai Gs", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="Achilles", dataType="rnai", fileType="genomicMatrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)


id <- "syn4598413"
fn <- "Achilles_shRNA_LINCS.csv"
achillesshRNAObj <- synGet(id)

achillesshRNA <- fread(getFileLocation(achillesshRNAObj), skip=2, data.table=FALSE,
                    select=c("Name", "Description", achillesLineNames)) %>%
  rename(MCF7=MCF7_BREAST, A549=A549_LUNG)

write.csv(achillesshRNA, file=fn)
f <- File(fn, name="Achilles QC rnai", parentId=proj@properties$id)
synSetAnnotations(f) <- list(source="Achilles", dataType="rnai", fileType="genomicMatrix")
act <- Activity(used=id, executed=thisScript)
generatedBy(f) <- act
f <- synStore(f)


