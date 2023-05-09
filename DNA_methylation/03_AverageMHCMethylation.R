# CALCLUATE AVERAGE MHCI METHYLATION:
# 10Oct22


# SET UP ENVIRONMENT & LOAD DATA:
# NB, forgot to actually save the README belonging to this dataset when I saved it in the previous script, so re-create in order to save

library(plyr)

load("AllEpicOrganoid_WorkingDataset_10Oct22.RData")
readme.epic <- data.frame(Dataset=c("samples","beta","mval"),Description=c("All QCed, filtered EPIC array samples (07Oct22)","Corresponding ComBat corrected beta values","Corresponding M values"))

readme.epic

gene.info <- read.csv("/home/fp215/rds/hpc-work/REFERENCE/EPIC_ARRAY/EPIC_ensembl_gene_annotation_RE.csv")[,c("IlmnID","CHR","MAPINFO","Gene.stable.ID","Gene.name","Gene.type","CpG_in")]


# MHCI GENES OF INTEREST:

goi <- c("NLRC5","HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","TAP1","TAP2","PSMB8","PSMB9","B2M","IRF1")

mhc <- gene.info[which(gene.info$Gene.name %in% goi),]

length(unique(mhc$IlmnID))


# EXTRACT ALL NLRC5/MHCI PROBES:

beta.mhc <- subset(beta, rownames(beta) %in% mhc$IlmnID)

nrow(beta.mhc)


# EXTRACT THE NLRC5 PROMPOTER PROBES:

nlrc5 <- subset(gene.info, Gene.name=="NLRC5" & CpG_in=="promoter")
beta.nlrc5 <- subset(beta, rownames(beta) %in% nlrc5$IlmnID)

ncol(beta.nlrc5)


# CALCULATE AVERAGE METHYLATION FOR ALL SAMPLES:

samples$Avg.MHCI.CpG <- colMeans(beta.mhc)


# RESAVE UPDATED DATASET:

update <- data.frame(Dataset=c("mhc","beta.mhc","beta.nlrc5","gene.info"),Description=c("MHCI genes of interest","Beta values for MHCI genes of interest","Beta values for NLRC5 promoter CpGs","CpG gene annotation info"))
readme.epic <- rbind(readme.epic,update)
save(samples,beta,mval,mhc,beta.mhc,beta.nlrc5,gene.info,readme.epic,file="AllEpicOrganoid_WorkingDataset_10Oct22.RData")


# SESSION INFO:

sessionInfo()


