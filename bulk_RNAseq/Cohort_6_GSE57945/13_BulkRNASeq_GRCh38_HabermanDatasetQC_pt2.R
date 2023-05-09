# QC OF HABERMAN GSE57945 RAW BULK RNASEQ COUNTS - Part 2
# Extract normalised counts (NB, no batch correction) & verify gender
# 17Jan22 (fp215)
# 
# Folders required:
#     - NORMALISED_COUNTS_QC
#       
# Files required:
#     - Raw counts dataset: /Users/fp215/Documents/DATASETS/RNASEQ/Haberman_GSE57945/GSE57945_COUNTS/RAW_COUNTS_QC/GSE57945_RawCounts.RData
#     - Raw counts DESeq2 dataset: /Users/fp215/Documents/DATASETS/RNASEQ/Haberman_GSE57945/GSE57945_COUNTS/RAW_COUNTS_QC/GSE57945_Raw_dds.RData


#######################
# WORKING ENVIRONMENT #
#######################

setwd("/Users/fp215/Documents/DATASETS/RNASEQ/Haberman_GSE57945")
library(DESeq2)
library(ClassDiscovery)
library(factoextra)
library(clustertend)
library(NbClust)
library(rafalib)
library(dendextend)
require(gridExtra)


#############
# LOAD DATA #
#############

load("GSE57945_RawCounts.RData")

readme.raw

load("GSE57945_Raw_dds.RData")

dds


#############################
# EXTRACT NORMALISED COUNTS #
#############################

norm <- counts(dds, normalized=T)

# Save normalised dataset:

readme.norm <- data.frame(Dataset=c("norm","pheno"),Description=c("Haberman GSE57945 DESeq2 normalised counts (n=322)","Haberman GSE57945 sample info (n=322)"))
save(norm,pheno,readme.norm,file="GSE57945_DESeq2NormalisedCounts.RData")


################
# CHECK GENDER #
################

# X: XIST (ENSG00000229807)
# Y: RPS4Y1 (ENSG00000129824)
# Y: EIF1AY (ENSG00000198692)
# Y: DDX3Y (ENSG00000067048)
# Y: KDM5D (ENSG00000012817)

# Subset the normalised counts to the sex-linked genes above:

sex.genes <- c("ENSG00000229807","ENSG00000129824","ENSG00000198692","ENSG00000067048","ENSG00000012817")
norm.sex <- subset(norm, rownames(norm) %in% sex.genes)

dim(norm.sex)

# Cluster plot of sex-linked genes:

Diagnosis <- c("skyblue","grey75","yellow")
Diagnosis <- Diagnosis[as.numeric(as.factor(pheno$Diagnosis))]
Sex <- c("darkred","darkblue")
Sex <- Sex[as.numeric(as.factor(pheno$Gender))]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

pdf("NORMALISED_COUNTS_QC/GSE57945_SexLinkedGenes_EuclideanAverageClustering.pdf", width=50, height=10)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=cbind(Diagnosis,Sex), dend=dend, rowLabels=c("Diagnosis","Sex"))
dev.off()


################
# SESSION INFO #
################

sessionInfo()
