# QC OF THE E-MTAB-5464 (HOWELL) RAW BULK RNASEQ COUNTS DATASET - Part 2
# Gender check, extract normalised counts & calculate average MHCI signature expression
# Data aligned to GRCh38
# 27Oct22 (fp215)
# 
# Folders required:
#     - NORMALISED_COUNTS
#     - NORMALISED_COUNTS_QC
#     - NORMALISED_COUNTS_QC/HOWELL_GRCh38
#
# Files required:
#     - Howell raw counts dataset: /home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/Howell_GRCh38_RawCounts.RData
#
# Functions used:
#     - Average expression function: /home/fp215/SCRIPTS/ADDENBROOKES_SCRIPTS/GENERAL_PURPOSE/AverageExpression_function.v1.R
#     - Plot themes & colour schemes: /home/fp215/SCRIPTS/ADDENBROOKES_SCRIPTS/GENERAL_PURPOSE/PlotTools.v3.4.R


#######################
# WORKING ENVIRONMENT #
#######################

setwd("/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER")
library(DESeq2)
library(sva)
library(ClassDiscovery)
library(factoextra)
library(clustertend)
library(NbClust)
library(rafalib)
library(dendextend)
require(gridExtra)
library(plyr)
source("/home/fp215/SCRIPTS/ADDENBROOKES_SCRIPTS/GENERAL_PURPOSE/AverageExpression_function.v1.R")
source("/home/fp215/SCRIPTS/ADDENBROOKES_SCRIPTS/GENERAL_PURPOSE/PlotTools.v3.4.R")

PCs_to_view <- 10


#################
# LOAD RAW DATA #
#################

# HOWELL DATASET:

load("RAW_COUNTS/Howell_GRCh38_RawCounts.RData")

readme.howell.raw


######################################
# RAW DATA: TISSUE SPECIFIC DATASETS #
######################################

# Also update the saved raw datasets to include these subsets

# HOWELL DATASET:
# TI, SC & AC samples

pheno.howell.ti <- subset(pheno.howell, tissue=="TI")
raw.howell.ti <- raw.howell[,pheno.howell.ti$rnaseq.id]
pheno.howell.sc <- subset(pheno.howell, tissue=="SC")
raw.howell.sc <- raw.howell[,pheno.howell.sc$rnaseq.id]
pheno.howell.ac <- subset(pheno.howell, tissue=="AC")
raw.howell.ac <- raw.howell[,pheno.howell.ac$rnaseq.id]

update <- data.frame(Dataset=c("raw.howell.ti","pheno.howell.ti","raw.howell.sc","pheno.howell.sc","raw.howell.ac","pheno.howell.ac"),Description=c("raw.howell subset to TI samples only (n=32)","TI sample info","raw.howell subset to SC samples only (n=32)","SC sample info","raw.howell subset to AC samples only (n=12)","AC samples only"))
readme.howell.raw <- rbind(readme.howell.raw,update)
save(raw.howell,pheno.howell,raw.howell.ti,pheno.howell.ti,raw.howell.sc,pheno.howell.sc,raw.howell.ac,pheno.howell.ac,info,readme.howell.raw,file="RAW_COUNTS/Howell_GRCh38_RawCounts.RData")


##############################
# RAW DATA: HEAT SCREE PLOTS #
##############################

# Subsetting each dataset to only the genes with a median of at least 5 counts across all samples

# HOWELL DATASETS:

pca <- raw.howell.ti[(apply(raw.howell.ti[,c(1:ncol(raw.howell.ti))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.howell.ti[c("sex","diagnosis")]
meta_continuous <- pheno.howell.ti[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("RAW_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_TI_RawCounts_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

pca <- raw.howell.sc[(apply(raw.howell.sc[,c(1:ncol(raw.howell.sc))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.howell.sc[c("sex","diagnosis")]
meta_continuous <- pheno.howell.sc[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("RAW_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_SC_RawCounts_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

pca <- raw.howell.ac[(apply(raw.howell.ac[,c(1:ncol(raw.howell.ac))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.howell.ac[c("sex","diagnosis")]
meta_continuous <- pheno.howell.ac[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("RAW_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_AC_RawCounts_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()


##########################
# CREATE DESEQ2 DATASETS #
##########################

# HOWELL DATASETS:

factors.howell.ti <- pheno.howell.ti[,c("rnaseq.id","case.no","diagnosis","sex","age")]
rownames(factors.howell.ti) <- factors.howell.ti$rnaseq.id
factors.howell.ti <- as.matrix(factors.howell.ti)

howell.dds.ti <- DESeqDataSetFromMatrix(countData=raw.howell.ti, colData=factors.howell.ti, design=~diagnosis)
howell.dds.ti$diagnosis <- relevel(howell.dds.ti$diagnosis, ref="Control")
howell.dds.ti <- howell.dds.ti[rowSums(counts(howell.dds.ti))>1,]
howell.vst.ti <- vst(howell.dds.ti,blind=F)
howell.dds.ti <- DESeq(howell.dds.ti)

factors.howell.sc <- pheno.howell.sc[,c("rnaseq.id","case.no","diagnosis","sex","age")]
rownames(factors.howell.sc) <- factors.howell.sc$rnaseq.id
factors.howell.sc <- as.matrix(factors.howell.sc)

howell.dds.sc <- DESeqDataSetFromMatrix(countData=raw.howell.sc, colData=factors.howell.sc, design=~diagnosis)
howell.dds.sc$diagnosis <- relevel(howell.dds.sc$diagnosis, ref="Control")
howell.dds.sc <- howell.dds.sc[rowSums(counts(howell.dds.sc))>1,]
howell.vst.sc <- vst(howell.dds.sc,blind=F)
howell.dds.sc <- DESeq(howell.dds.sc)

factors.howell.ac <- pheno.howell.ac[,c("rnaseq.id","case.no","diagnosis","sex","age")]
rownames(factors.howell.ac) <- factors.howell.ac$rnaseq.id
factors.howell.ac <- as.matrix(factors.howell.ac)

howell.dds.ac <- DESeqDataSetFromMatrix(countData=raw.howell.ac, colData=factors.howell.ac, design=~diagnosis)
howell.dds.ac$diagnosis <- relevel(howell.dds.ac$diagnosis, ref="Control")
howell.dds.ac <- howell.dds.ac[rowSums(counts(howell.dds.ac))>1,]
howell.vst.ac <- vst(howell.dds.ac,blind=F)
howell.dds.ac <- DESeq(howell.dds.ac)

save(howell.dds.ti,howell.vst.ti,factors.howell.ti,howell.dds.sc,howell.vst.sc,factors.howell.sc,howell.dds.ac,howell.vst.ac,factors.howell.ac,file="RAW_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_ByTissue_RawDataQCInput.RData")


#############################
# EXTRACT NORMALISED COUNTS #
#############################

# HOWELL DATASETS:

howell.ti <- counts(howell.dds.ti, normalized=T)
howell.sc <- counts(howell.dds.sc, normalized=T)
howell.ac <- counts(howell.dds.ac, normalized=T)

readme.howell.norm <- data.frame(Dataset=c("howell.ti","pheno.howell.ti","howell.sc","pheno.howell.sc","howell.ac","pheno.howell.ac","info"),Description=c("E-MTAB-5464 (Howell) normalised TI counts (n=32)","Howell TI sample info, updated 29Sep22","E-MTAB-5464 (Howell) normalised SC counts (n=32)","Howell SC sample info, updated 29Sep22","E-MTAB-5464 (Howell) normalised AC counts (n=12)","Howell AC sample info, updated 29Sep22","Ensembl GRCh38 release 99 gene annotation info"))
save(howell.ti,pheno.howell.ti,howell.sc,pheno.howell.sc,howell.ac,pheno.howell.ac,info,readme.howell.norm,file="NORMALISED_COUNTS/Howell_GRCh38_ByTissue_NormalisedCounts.RData")


################################
# NLRC5/MHCI GENES OF INTEREST #
################################

goi <- data.frame(Gene=c("NLRC5","HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","TAP1","TAP2","PSMB8","PSMB9","B2M","IRF1"),Ensembl.ID=c("ENSG00000140853","ENSG00000206503","ENSG00000234745","ENSG00000204525","ENSG00000204592","ENSG00000204642","ENSG00000204632","ENSG00000168394","ENSG00000204267","ENSG00000204264","ENSG00000240065","ENSG00000166710","ENSG00000125347"))

goi


#################################################
# NORMALISED DATA: CALCULATE AVERAGE EXPRESSION #
#################################################

# Input fields:
#	- assay.data : matrix of normalised counts
#	- features : list of vectors of features (each vector containing a gene set of interest)
#	- nbin : number of bins of aggregate expression levels required for all analysed features (default = 50)
#	- controls : number of control features to be selected from the same bin per analysed feature (default = 100)
#	- name : name preceding each average expression score produced for each vector of features (unimportant)
#	- seed : random seed (default = 1)
#	- genekey : a dataframe containing a gene ID : gene symbol key (default = NULL)
#	- info : any additional info to be included in the output object (default = NULL)


# HOWELL DATASETS:

howell.ti.avexprs <- AverageExpression(howell.ti, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.howell.ti)
howell.sc.avexprs <- AverageExpression(howell.sc, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.howell.sc)
howell.ac.avexprs <- AverageExpression(howell.ac, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.howell.ac)

# UPDATE & SAVE METADATA FOR ALL DATASETS:

ti.avg <- data.frame(rnaseq.id=rownames(howell.ti.avexprs$AverageExpression),Avg.MHCI.exprs=(howell.ti.avexprs$AverageExpression)$Avg.1)
pheno.howell.ti <- join(pheno.howell.ti, ti.avg, type="left", match="all")
sc.avg <- data.frame(rnaseq.id=rownames(howell.sc.avexprs$AverageExpression),Avg.MHCI.exprs=(howell.sc.avexprs$AverageExpression)$Avg.1)
pheno.howell.sc <- join(pheno.howell.sc, sc.avg, type="left", match="all")
ac.avg <- data.frame(rnaseq.id=rownames(howell.ac.avexprs$AverageExpression),Avg.MHCI.exprs=(howell.ac.avexprs$AverageExpression)$Avg.1)
pheno.howell.ac <- join(pheno.howell.ac, ac.avg, type="left", match="all")
save(howell.ti,pheno.howell.ti,howell.sc,pheno.howell.sc,howell.ac,pheno.howell.ac,info,readme.howell.norm,file="NORMALISED_COUNTS/Howell_GRCh38_ByTissue_NormalisedCounts.RData")


##############################
# NORMALISED DATA: PCA PLOTS #
##############################

# HOWELL DATASETS:

pca <- prcomp(t(howell.ti))
pca <- as.data.frame(pca$x)[,1:10]
pca$rnaseq.id <- rownames(pca)
pca <- join(pheno.howell.ti, pca, type="left", match="all")
write.table(pca, "NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_TI_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_TI_NormalisedCounts_PC1vPC2.pdf", width=10, height=5)
p1 <- ggplot(pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis6,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=2)
dev.off()

pca <- prcomp(t(howell.sc))
pca <- as.data.frame(pca$x)[,1:10]
pca$rnaseq.id <- rownames(pca)
pca <- join(pheno.howell.sc, pca, type="left", match="all")
write.table(pca, "NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_SC_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_SC_NormalisedCounts_PC1vPC2.pdf", width=10, height=5)
p1 <- ggplot(pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis6,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=2)
dev.off()

pca <- prcomp(t(howell.ac))
pca <- as.data.frame(pca$x)[,1:10]
pca$rnaseq.id <- rownames(pca)
pca <- join(pheno.howell.ac, pca, type="left", match="all")
write.table(pca, "NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_AC_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_AC_NormalisedCounts_PC1vPC2.pdf", width=10, height=5)
p1 <- ggplot(pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis6,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=2)
dev.off()


#####################################
# NORMALISED DATA: HEAT SCREE PLOTS #
#####################################

# Subsetting each dataset to only the genes with a median of at least 5 counts across all samples

# HOWELL DATASETS:

pca <- howell.ti[(apply(howell.ti[,c(1:ncol(howell.ti))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.howell.ti[c("sex","diagnosis")]
meta_continuous <- pheno.howell.ti[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_TI_NormalisedCounts_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

pca <- howell.sc[(apply(howell.sc[,c(1:ncol(howell.sc))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.howell.sc[c("sex","diagnosis")]
meta_continuous <- pheno.howell.sc[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_SC_NormalisedCounts_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

pca <- howell.ac[(apply(howell.ac[,c(1:ncol(howell.ac))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.howell.ac[c("sex","diagnosis")]
meta_continuous <- pheno.howell.ac[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_AC_NormalisedCounts_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()


################
# CHECK GENDER #
################

# X: XIST (ENSG00000229807)
# Y: RPS4Y1 (ENSG00000129824)
# Y: EIF1AY (ENSG00000198692)
# Y: DDX3Y (ENSG00000067048)
# Y: KDM5D (ENSG00000012817)

sex.genes <- c("ENSG00000229807","ENSG00000129824","ENSG00000198692","ENSG00000067048","ENSG00000012817")

# HOWELL DATASETS:

norm.sex <- subset(howell.ti, rownames(howell.ti) %in% sex.genes)

Diagnosis <- c("skyblue","grey75")
Diagnosis <- Diagnosis[as.numeric(as.factor(pheno.howell.ti$diagnosis))]
Sex <- c("darkred","darkblue")
Sex <- Sex[as.numeric(as.factor(pheno.howell.ti$sex))]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

pdf("NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_TI_SexLinkedGenes_EuclideanAverageClustering.pdf", width=20, height=11)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=cbind(Diagnosis,Sex), dend=dend, rowLabels=c("Diagnosis","Sex"))
dev.off()

norm.sex <- subset(howell.sc, rownames(howell.sc) %in% sex.genes)

Diagnosis <- c("skyblue","grey75")
Diagnosis <- Diagnosis[as.numeric(as.factor(pheno.howell.sc$diagnosis))]
Sex <- c("darkred","darkblue")
Sex <- Sex[as.numeric(as.factor(pheno.howell.sc$sex))]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

pdf("NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_SC_SexLinkedGenes_EuclideanAverageClustering.pdf", width=20, height=11)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=cbind(Diagnosis,Sex), dend=dend, rowLabels=c("Diagnosis","Sex"))
dev.off()

norm.sex <- subset(howell.ac, rownames(howell.ac) %in% sex.genes)

Diagnosis <- c("skyblue","grey75")
Diagnosis <- Diagnosis[as.numeric(as.factor(pheno.howell.ac$diagnosis))]
Sex <- c("darkred","darkblue")
Sex <- Sex[as.numeric(as.factor(pheno.howell.ac$sex))]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

pdf("NORMALISED_COUNTS_QC/HOWELL_GRCh38/Howell_GRCh38_AC_SexLinkedGenes_EuclideanAverageClustering.pdf", width=15, height=11)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=cbind(Diagnosis,Sex), dend=dend, rowLabels=c("Diagnosis","Sex"))
dev.off()


################
# SESSION INFO #
################

sessionInfo()
rm(list=ls())

