# QC OF THE NGS-K.NAYAK-40706_RE (E-MTAB-11548) RAW BULK RNASEQ COUNTS DATASET - Part 2
# Removal of differentiated samples before normalisation of the RE dataset
# Gender check, extract normalised counts & calculate average MHCI signature expression
# Data aligned to GRCh38
# 27Oct22 (fp215)
# 
# Folders required:
#     - NORMALISED_COUNTS
#     - NORMALISED_COUNTS_QC
#     - NORMALISED_COUNTS_QC/RE_GRCh38
#
# Files required:
#     - RE raw counts dataset: /home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/RE_GRCh38_RawCounts.RData
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

# RE DATASET:

load("RAW_COUNTS/RE_GRCh38_RawCounts.RData")

readme.re.raw


######################################
# RAW DATA: TISSUE SPECIFIC DATASETS #
######################################

# Also update the saved raw datasets to include these subsets

# RE DATASET:
# TI & SC samples
# NB, dataset also includes TNFa/IFNg treated samples

pheno.re.ti <- subset(pheno.re, tissue=="TI")
raw.re.ti <- raw.re[,pheno.re.ti$rnaseq.id]
pheno.re.sc <- subset(pheno.re, tissue=="SC")
raw.re.sc <- raw.re[,pheno.re.sc$rnaseq.id]

update <- data.frame(Dataset=c("raw.re.ti","pheno.re.ti","raw.re.sc","pheno.re.sc"),Description=c("raw.re subset to TI samples only (n=24)","TI sample info","raw.re subset to SC samples only (n=6)","SC sample info"))
readme.re.raw <- rbind(readme.re.raw,update)
save(raw.re.all,pheno.re.all,raw.re,pheno.re,info,raw.re.ti,pheno.re.ti,raw.re.sc,pheno.re.sc,readme.re.raw,file="RAW_COUNTS/RE_GRCh38_RawCounts.RData")


##############################
# RAW DATA: HEAT SCREE PLOTS #
##############################

# Subsetting each dataset to only the genes with a median of at least 5 counts across all samples

# RE UNDIFFERENTIATED DATASETS:
# NB, there are only 6 SC samples, so it doesn't make sense to create a heat scree plot for these

pca <- raw.re.ti[(apply(raw.re.ti[,c(1:ncol(raw.re.ti))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.re.ti[c("sex","diagnosis","treatment")]
meta_continuous <- pheno.re.ti[c("age","passage")]
meta_continuous$passage <- as.numeric(gsub("P","",meta_continuous$passage))
colnames(meta_categorical) <- c("Sex","Diagnosis","Treatment")
colnames(meta_continuous) <- c("Age","Passage")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_UndifferentiatedTI_RawCounts_HeatScree.pdf", width=10, height=7)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()


##########################
# CREATE DESEQ2 DATASETS #
##########################

# RE UNDIFFERENTIATED DATASETS:

re.factors <- pheno.re[,c("rnaseq.id","case.no","diagnosis","sex","age","passage","tissue","treatment")]
re.factors$group <- paste0(re.factors$tissue,".",re.factors$diagnosis)
re.factors$group2 <- paste0(re.factors$group,".",re.factors$treatment)
re.factors$group3 <- paste0(re.factors$diagnosis,".",re.factors$treatment)
rownames(re.factors) <- re.factors$rnaseq.id
re.factors$diagnosis <- as.factor(re.factors$diagnosis)

re.dds <- DESeqDataSetFromMatrix(countData=raw.re, colData=re.factors, design=~group2)
re.dds$diagnosis <- relevel(re.dds$diagnosis, ref="Control")
re.dds <- re.dds[rowSums(counts(re.dds))>1,]
re.vst <- vst(re.dds,blind=F)
re.dds <- DESeq(re.dds)

re.factors.ti <- pheno.re.ti[,c("rnaseq.id","case.no","diagnosis","sex","age","passage","treatment")]
re.factors.ti$group <- paste0(re.factors.ti$diagnosis,".",re.factors.ti$treatment)
rownames(re.factors.ti) <- re.factors.ti$rnaseq.id
re.factors.ti$diagnosis <- as.factor(re.factors.ti$diagnosis)

re.dds.ti <- DESeqDataSetFromMatrix(countData=raw.re.ti, colData=re.factors.ti, design=~group)
re.dds.ti$diagnosis <- relevel(re.dds.ti$diagnosis, ref="Control")
re.dds.ti <- re.dds.ti[rowSums(counts(re.dds.ti))>1,]
re.vst.ti <- vst(re.dds.ti,blind=F)
re.dds.ti <- DESeq(re.dds.ti)

re.factors.sc <- pheno.re.sc[,c("rnaseq.id","case.no","diagnosis","sex","age","passage","treatment")]
re.factors.sc$group <- paste0(re.factors.sc$diagnosis,".",re.factors.sc$treatment)
rownames(re.factors.sc) <- re.factors.sc$rnaseq.id
re.factors.sc$diagnosis <- as.factor(re.factors.sc$diagnosis)

re.dds.sc <- DESeqDataSetFromMatrix(countData=raw.re.sc, colData=re.factors.sc, design=~1)
re.dds.sc <- re.dds.sc[rowSums(counts(re.dds.sc))>1,]
re.vst.sc <- vst(re.dds.sc,blind=F)
re.dds.sc <- DESeq(re.dds.sc)

save(re.dds.ti,re.vst.ti,re.factors.ti,re.dds.sc,re.vst.sc,re.factors.sc,re.dds,re.vst,re.factors,file="RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_ByTissue_RawDataQCInput.RData")


#############################
# EXTRACT NORMALISED COUNTS #
#############################

# RE UNDIFFERENTIATED DATASETS:

re <- counts(re.dds, normalized=T)
re.ti <- counts(re.dds.ti, normalized=T)
re.sc <- counts(re.dds.sc, normalized=T)

readme.re.norm <- data.frame(Dataset=c("re","pheno.re","re.ti","pheno.re.ti","re.sc","pheno.re.sc","info"),Description=c("NGS-K.Nayak-40706 (RE) undifferentiated normalised counts (n=30)","RE sample info, updated 24Oct22","NGS-K.Nayak-40706 (RE) undifferentiated normalised TI counts (n=24)","RE TI sample info, updated 24Oct22","NGS-K.Nayak-40706 (RE) undifferentiated normalised SC counts (n=6)","RE SC sample info, updated 24Oct22","Ensembl GRCh38 release 99 gene annotation info"))
save(re,pheno.re,re.ti,pheno.re.ti,re.sc,pheno.re.sc,info,readme.re.norm,file="NORMALISED_COUNTS/RE_GRCh38_ByTissue_NormalisedCounts.RData")


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


# RE UNDIFFERENTIATED DATASETS:

re.avexprs <- AverageExpression(re, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.re)
re.ti.avexprs <- AverageExpression(re.ti, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.re.ti)
re.sc.avexprs <- AverageExpression(re.sc, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.re.sc)

# UPDATE & SAVE METADATA FOR ALL DATASETS:

re.avg <- data.frame(rnaseq.id=rownames(re.avexprs$AverageExpression),Avg.MHCI.exprs=(re.avexprs$AverageExpression)$Avg.1)
pheno.re <- join(pheno.re, re.avg, type="left", match="all")
ti.avg <- data.frame(rnaseq.id=rownames(re.ti.avexprs$AverageExpression),Avg.MHCI.exprs=(re.ti.avexprs$AverageExpression)$Avg.1)
pheno.re.ti <- join(pheno.re.ti, ti.avg, type="left", match="all")
sc.avg <- data.frame(rnaseq.id=rownames(re.sc.avexprs$AverageExpression),Avg.MHCI.exprs=(re.sc.avexprs$AverageExpression)$Avg.1)
pheno.re.sc <- join(pheno.re.sc, sc.avg, type="left", match="all")
save(re,pheno.re,re.ti,pheno.re.ti,re.sc,pheno.re.sc,info,readme.re.norm,file="NORMALISED_COUNTS/RE_GRCh38_ByTissue_NormalisedCounts.RData")
write.table(pheno.re, "NORMALISED_COUNTS/RE_GRCh38_UndifferentiatedSamples_AverageMHCIExpression.txt", row.names=F, sep="\t", quote=F)


##############################
# NORMALISED DATA: PCA PLOTS #
##############################

# RE UNDIFFERENTIATED DATASETS:

r.pca <- prcomp(t(re))
r.pca <- as.data.frame(r.pca$x)[,1:10]
r.pca$rnaseq.id <- rownames(r.pca)
r.pca <- join(pheno.re, r.pca, type="left", match="all")
write.table(r.pca, "NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_UndifferentiatedALL_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_UndifferentiatedALL_NormalisedCounts_PC1vPC2.pdf", width=15, height=5)
p1 <- ggplot(r.pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(r.pca, aes(PC1, PC2, fill=tissue))
p3 <- ggplot(r.pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite2,p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=3)
dev.off()

pca <- prcomp(t(re.ti))
pca <- as.data.frame(pca$x)[,1:10]
pca$rnaseq.id <- rownames(pca)
pca <- join(pheno.re.ti, pca, type="left", match="all")
write.table(pca, "NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_UndifferentiatedTI_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_UndifferentiatedTI_NormalisedCounts_PC1vPC2.pdf", width=10, height=10)
p1 <- ggplot(pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(pca, aes(PC1, PC2, fill=treatment))
p4 <- ggplot(pca, aes(PC1, PC2, fill=passage))
p6 <- ggplot(pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_treat,p4+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Passage",values=c(P2="red",P3="darkorange",P4="gold",P8="green4",P9="mediumblue",P12="magenta4")),p6+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=3)
dev.off()

pca <- prcomp(t(re.sc))
pca <- as.data.frame(pca$x)
pca$rnaseq.id <- rownames(pca)
pca <- join(pheno.re.sc, pca, type="left", match="all")
write.table(pca, "NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_UndifferentiatedSC_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_UndifferentiatedSC_NormalisedCounts_PC1vPC2.pdf", width=10, height=10)
p1 <- ggplot(pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(pca, aes(PC1, PC2, fill=treatment))
p4 <- ggplot(pca, aes(PC1, PC2, fill=passage))
p6 <- ggplot(pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_treat,p4+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Passage",values=c(P2="red",P3="darkorange",P4="gold",P9="mediumblue",P12="magenta4")),p6+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=3)
dev.off()


#####################################
# NORMALISED DATA: HEAT SCREE PLOTS #
#####################################

# Subsetting each dataset to only the genes with a median of at least 5 counts across all samples

# RE UNDIFFERENTIATED DATASETS:
# NB, there are only 6 SC samples, so it doesn't make sense to create a heat scree plot for these

pca <- re[(apply(re[,c(1:ncol(re))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.re[c("sex","diagnosis","treatment","tissue")]
meta_continuous <- pheno.re[c("age","passage")]
meta_continuous$passage <- as.numeric(gsub("P","",meta_continuous$passage))
colnames(meta_categorical) <- c("Sex","Diagnosis","Treatment","Tissue")
colnames(meta_continuous) <- c("Age","Passage")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_UndifferentiatedALL_NormalisedCounts_HeatScree.pdf", width=10, height=7)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

pca <- re.ti[(apply(re.ti[,c(1:ncol(re.ti))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.re.ti[c("sex","diagnosis","treatment")]
meta_continuous <- pheno.re.ti[c("age","passage")]
meta_continuous$passage <- as.numeric(gsub("P","",meta_continuous$passage))
colnames(meta_categorical) <- c("Sex","Diagnosis","Treatment")
colnames(meta_continuous) <- c("Age","Passage")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_UndifferentiatedTI_NormalisedCounts_HeatScree.pdf", width=10, height=7)
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

# RE UNDIFFERENTIATED DATASETS:

norm.sex <- subset(re.ti, rownames(re.ti) %in% sex.genes)

Diagnosis <- c("skyblue","grey75")
Diagnosis <- Diagnosis[as.numeric(as.factor(pheno.re.ti$diagnosis))]
Sex <- c("darkred","darkblue")
Sex <- Sex[as.numeric(as.factor(pheno.re.ti$sex))]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

pdf("NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_UndifferentiatedTI_SexLinkedGenes_EuclideanAverageClustering.pdf", width=20, height=11)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=cbind(Diagnosis,Sex), dend=dend, rowLabels=c("Diagnosis","Sex"))
dev.off()

norm.sex <- subset(re.sc, rownames(re.sc) %in% sex.genes)

Diagnosis <- c("skyblue","grey75")
Diagnosis <- Diagnosis[as.numeric(as.factor(pheno.re.sc$diagnosis))]
Sex <- c("darkred","darkblue")
Sex <- Sex[as.numeric(as.factor(pheno.re.sc$sex))]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

pdf("NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_UndifferentiatedSC_SexLinkedGenes_EuclideanAverageClustering.pdf", width=10, height=11)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=cbind(Diagnosis,Sex), dend=dend, rowLabels=c("Diagnosis","Sex"))
dev.off()


################
# SESSION INFO #
################

sessionInfo()
rm(list=ls())

