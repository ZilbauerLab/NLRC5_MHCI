# LOAD UP DATASET, PULL OUT M VALUES & RE-SAVE
# 10Oct22


# SET ENVIRONMENT & LOAD DATA:

library(plyr)

load("AllEPICOrganoid_AnalysisReadyDataset_10Oct22.RData")

ls()

readme.epic

convert <- function(x){log2(x/(1-x))}


# CONVERT BETA VALUES TO M VALUES (WHOLE DATASET):

all.mval <- structure(sapply(all.beta, convert), dim=dim(all.beta))
rownames(all.mval) <- rownames(all.beta)
colnames(all.mval) <- colnames(all.beta)

dim(all.mval)


# FINAL DATASET:


beta <- all.beta[,samples$array.id]
mval <- all.mval[,samples$array.id]


# RE-SAVE FULL DATASET & CREATE AN ANALYSIS-READY DATASET:

update <- data.frame(Dataset=c("all.mval","beta","mval"),Description=c("All RE QCed EPIC methylation array M values (10Oct22)","Final analysis ready filtered beta values (corresponding to samples, 10Oct22)","Final analysis ready filtered M values (corresponding to samples, 10Oct22)"))
readme.epic <- rbind(readme.epic,update)
save(epic.samples,all.beta,re.meta,all.ok.samples,samples,beta,mval,all.mval,readme.epic,file="AllEPICOrganoid_AnalysisReadyDataset_10Oct22.RData")

readme.epic <- data.frame(Dataset=c("samples","beta","mval"),Description=c("All QCed, filtered EPIC array samples (07Oct22)","Corresponding ComBat corrected beta values","Corresponding M values"))
save(samples,beta,mval,file="AllEpicOrganoid_WorkingDataset_10Oct22.RData")


# Session Info:

sessionInfo()





