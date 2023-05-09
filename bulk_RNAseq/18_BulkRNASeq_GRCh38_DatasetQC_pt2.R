# QC OF THE E-MTAB-5464 (HOWELL), NGS-K.NAYAK-40706_RE (RE) & NGS-K.NAYAK-40750_AF (AF) RAW BULK RNASEQ COUNTS DATASETS - Part 2
# Removal of differentiated samples before normalisation of the RE dataset
# Gender check, extract normalised counts & calculate average MHCI signature expression
# Data aligned to GRCh38
# 05Sep22 (fp215)
# 27Oct22 (fp215) - updates using separated AF & RE datasets + updated sample info. Differentiated samples already removed from the RE dataset.
# 
# Folders required:
#     - NORMALISED_COUNTS
#     - NORMALISED_COUNTS_QC
#     - NORMALISED_COUNTS_QC/HOWELL_GRCh38
#     - NORMALISED_COUNTS_QC/RE_GRCh38
#     - NORMALISED_COUNTS_QC/AF_GRCh38
#
# Files required:
#     - Howell raw counts dataset: /home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/Howell_GRCh38_RawCounts.RData
#     - RE raw counts dataset: /home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/RE_GRCh38_RawCounts.RData
#     - AF raw counts dataset: /home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/AF_GRCh38_RawCounts.RData
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

# RE DATASET:

load("RAW_COUNTS/RE_GRCh38_RawCounts.RData")

readme.re.raw

# AF DATASET:

load("RAW_COUNTS/AF_GRCh38_RawCounts.RData")

readme.af.raw


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

# AF DATASET:
# TI, SC & DUO samples

pheno.af.ti <- subset(pheno.af, tissue=="TI")
raw.af.ti <- raw.af[,pheno.af.ti$rnaseq.id]
pheno.af.sc <- subset(pheno.af, tissue=="SC")
raw.af.sc <- raw.af[,pheno.af.sc$rnaseq.id]
pheno.af.duo <- subset(pheno.af, tissue=="DUO")
raw.af.duo <- raw.af[,pheno.af.duo$rnaseq.id]

update <- data.frame(Dataset=c("raw.af.ti","pheno.af.ti","raw.af.sc","pheno.af.sc","raw.af.duo","pheno.af.duo"),Description=c("raw.af subset to TI samples only (n=30)","TI sample info","raw.af subset to SC samples only (n=30)","SC sample info","raw.af subset to DUO samples only (n=30)","DUO sample info"))
readme.af.raw <- rbind(readme.af.raw,update)
save(raw.af,pheno.af,info,raw.af.ti,pheno.af.ti,raw.af.sc,pheno.af.sc,raw.af.duo,pheno.af.duo,readme.af.raw,file="RAW_COUNTS/AF_GRCh38_RawCounts.RData")


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


# AF DATASETS:
# NB, all samples are from passage 3

pca <- raw.af.ti[(apply(raw.af.ti[,c(1:ncol(raw.af.ti))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.af.ti[c("sex","diagnosis","treatment")]
meta_continuous <- pheno.af.ti[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis","Treatment")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("RAW_COUNTS_QC/AF_GRCh38/AF_GRCh38_TI_RawCounts_HeatScree.pdf", width=10, height=7)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

pca <- raw.af.sc[(apply(raw.af.sc[,c(1:ncol(raw.af.sc))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.af.sc[c("sex","diagnosis","treatment")]
meta_continuous <- pheno.af.sc[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis","Treatment")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("RAW_COUNTS_QC/AF_GRCh38/AF_GRCh38_SC_RawCounts_HeatScree.pdf", width=10, height=7)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

pca <- raw.af.duo[(apply(raw.af.duo[,c(1:ncol(raw.af.duo))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.af.duo[c("sex","diagnosis","treatment")]
meta_continuous <- pheno.af.duo[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis","Treatment")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("RAW_COUNTS_QC/AF_GRCh38/AF_GRCh38_DUO_RawCounts_HeatScree.pdf", width=10, height=7)
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


# AF DATASETS:

af.factors <- pheno.af[,c("rnaseq.id","case.no","diagnosis","sex","age","passage","tissue","treatment")]
af.factors$group <- paste0(af.factors$tissue,".",af.factors$diagnosis)
af.factors$group2 <- paste0(af.factors$group,".",af.factors$treatment)
af.factors$group3 <- paste0(af.factors$diagnosis,".",af.factors$treatment)
rownames(af.factors) <- af.factors$rnaseq.id
af.factors$diagnosis <- as.factor(af.factors$diagnosis)

af.dds <- DESeqDataSetFromMatrix(countData=raw.af, colData=af.factors, design=~group2)
af.dds$diagnosis <- relevel(af.dds$diagnosis, ref="Control")
af.dds <- af.dds[rowSums(counts(af.dds))>1,]
af.vst <- vst(af.dds,blind=F)
af.dds <- DESeq(af.dds)

af.factors.ti <- pheno.af.ti[,c("rnaseq.id","case.no","diagnosis","sex","age","passage","treatment")]
af.factors.ti$group <- paste0(af.factors.ti$diagnosis,".",af.factors.ti$treatment)
rownames(af.factors.ti) <- af.factors.ti$rnaseq.id
af.factors.ti$diagnosis <- as.factor(af.factors.ti$diagnosis)

af.dds.ti <- DESeqDataSetFromMatrix(countData=raw.af.ti, colData=af.factors.ti, design=~group)
af.dds.ti$diagnosis <- relevel(af.dds.ti$diagnosis, ref="Control")
af.dds.ti <- af.dds.ti[rowSums(counts(af.dds.ti))>1,]
af.vst.ti <- vst(af.dds.ti,blind=F)
af.dds.ti <- DESeq(af.dds.ti)

af.factors.sc <- pheno.af.sc[,c("rnaseq.id","case.no","diagnosis","sex","age","passage","treatment")]
af.factors.sc$group <- paste0(af.factors.sc$diagnosis,".",af.factors.sc$treatment)
rownames(af.factors.sc) <- af.factors.sc$rnaseq.id
af.factors.sc$diagnosis <- as.factor(af.factors.sc$diagnosis)

af.dds.sc <- DESeqDataSetFromMatrix(countData=raw.af.sc, colData=af.factors.sc, design=~group)
af.dds.sc$diagnosis <- relevel(af.dds.sc$diagnosis, ref="Control")
af.dds.sc <- af.dds.sc[rowSums(counts(af.dds.sc))>1,]
af.vst.sc <- vst(af.dds.sc,blind=F)
af.dds.sc <- DESeq(af.dds.sc)

af.factors.duo <- pheno.af.duo[,c("rnaseq.id","case.no","diagnosis","sex","age","passage","treatment")]
af.factors.duo$group <- paste0(af.factors.duo$diagnosis,".",af.factors.duo$treatment)
rownames(af.factors.duo) <- af.factors.duo$rnaseq.id
af.factors.duo$diagnosis <- as.factor(af.factors.duo$diagnosis)

af.dds.duo <- DESeqDataSetFromMatrix(countData=raw.af.duo, colData=af.factors.duo, design=~group)
af.dds.duo$diagnosis <- relevel(af.dds.duo$diagnosis, ref="Control")
af.dds.duo <- af.dds.duo[rowSums(counts(af.dds.duo))>1,]
af.vst.duo <- vst(af.dds.duo,blind=F)
af.dds.duo <- DESeq(af.dds.duo)

save(af.dds.ti,af.vst.ti,af.factors.ti,af.dds.sc,af.vst.sc,af.factors.sc,af.dds.duo,af.vst.duo,af.factors.duo,af.dds,af.vst,af.factors,file="RAW_COUNTS_QC/AF_GRCh38/AF_GRCh38_ByTissue_RawDataQCInput.RData")


#############################
# EXTRACT NORMALISED COUNTS #
#############################

# HOWELL DATASETS:

howell.ti <- counts(howell.dds.ti, normalized=T)
howell.sc <- counts(howell.dds.sc, normalized=T)
howell.ac <- counts(howell.dds.ac, normalized=T)

readme.howell.norm <- data.frame(Dataset=c("howell.ti","pheno.howell.ti","howell.sc","pheno.howell.sc","howell.ac","pheno.howell.ac","info"),Description=c("E-MTAB-5464 (Howell) normalised TI counts (n=32)","Howell TI sample info, updated 29Sep22","E-MTAB-5464 (Howell) normalised SC counts (n=32)","Howell SC sample info, updated 29Sep22","E-MTAB-5464 (Howell) normalised AC counts (n=12)","Howell AC sample info, updated 29Sep22","Ensembl GRCh38 release 99 gene annotation info"))
save(howell.ti,pheno.howell.ti,howell.sc,pheno.howell.sc,howell.ac,pheno.howell.ac,info,readme.howell.norm,file="NORMALISED_COUNTS/Howell_GRCh38_ByTissue_NormalisedCounts.RData")


# RE UNDIFFERENTIATED DATASETS:

re <- counts(re.dds, normalized=T)
re.ti <- counts(re.dds.ti, normalized=T)
re.sc <- counts(re.dds.sc, normalized=T)

readme.re.norm <- data.frame(Dataset=c("re","pheno.re","re.ti","pheno.re.ti","re.sc","pheno.re.sc","info"),Description=c("NGS-K.Nayak-40706 (RE) undifferentiated normalised counts (n=30)","RE sample info, updated 24Oct22","NGS-K.Nayak-40706 (RE) undifferentiated normalised TI counts (n=24)","RE TI sample info, updated 24Oct22","NGS-K.Nayak-40706 (RE) undifferentiated normalised SC counts (n=6)","RE SC sample info, updated 24Oct22","Ensembl GRCh38 release 99 gene annotation info"))
save(re,pheno.re,re.ti,pheno.re.ti,re.sc,pheno.re.sc,info,readme.re.norm,file="NORMALISED_COUNTS/RE_GRCh38_ByTissue_NormalisedCounts.RData")


# AF DATASETS:

af <- counts(af.dds, normalized=T)
af.ti <- counts(af.dds.ti, normalized=T)
af.sc <- counts(af.dds.sc, normalized=T)
af.duo <- counts(af.dds.duo, normalized=T)

readme.af.norm <- data.frame(Dataset=c("af","pheno.af","af.ti","pheno.af.ti","af.sc","pheno.af.sc","af.duo","pheno.af.duo","info"),Description=c("NGS-K.Nayak-40750 (AF) normalised counts (n=90)","AF sample info, updated 24Oct22","NGS-K.Nayak-40750 (AF) normalised TI counts (n=30)","AF TI sample info, updated 24Oct22","NGS-K.Nayak-40750 (AF) normalised SC counts (n=30)","AF SC sample info, updated 24Oct22","NGS-K.Nayak-40750 (AF) normalised DUO counts (n=30)","AF DUO sample info, updated 24Oct22","Ensembl GRCh38 release 99 gene annotation info"))
save(af,pheno.af,af.ti,pheno.af.ti,af.sc,pheno.af.sc,af.duo,pheno.af.duo,info,readme.af.norm,file="NORMALISED_COUNTS/AF_GRCh38_ByTissue_NormalisedCounts.RData")


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


# RE UNDIFFERENTIATED DATASETS:

re.avexprs <- AverageExpression(re, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.re)
re.ti.avexprs <- AverageExpression(re.ti, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.re.ti)
re.sc.avexprs <- AverageExpression(re.sc, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.re.sc)


# AF DATASETS:

af.avexprs <- AverageExpression(af, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.af)
af.ti.avexprs <- AverageExpression(af.ti, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.af.ti)
af.sc.avexprs <- AverageExpression(af.sc, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.af.sc)
af.duo.avexprs <- AverageExpression(af.duo, list(goi$Ensembl.ID), nbin=50, controls=100, name="Avg.", seed=1, genekey=goi, info=pheno.af.duo)


# UPDATE & SAVE METADATA FOR ALL DATASETS:

ti.avg <- data.frame(rnaseq.id=rownames(howell.ti.avexprs$AverageExpression),Avg.MHCI.exprs=(howell.ti.avexprs$AverageExpression)$Avg.1)
pheno.howell.ti <- join(pheno.howell.ti, ti.avg, type="left", match="all")
sc.avg <- data.frame(rnaseq.id=rownames(howell.sc.avexprs$AverageExpression),Avg.MHCI.exprs=(howell.sc.avexprs$AverageExpression)$Avg.1)
pheno.howell.sc <- join(pheno.howell.sc, sc.avg, type="left", match="all")
ac.avg <- data.frame(rnaseq.id=rownames(howell.ac.avexprs$AverageExpression),Avg.MHCI.exprs=(howell.ac.avexprs$AverageExpression)$Avg.1)
pheno.howell.ac <- join(pheno.howell.ac, ac.avg, type="left", match="all")
save(howell.ti,pheno.howell.ti,howell.sc,pheno.howell.sc,howell.ac,pheno.howell.ac,info,readme.howell.norm,file="NORMALISED_COUNTS/Howell_GRCh38_ByTissue_NormalisedCounts.RData")

re.avg <- data.frame(rnaseq.id=rownames(re.avexprs$AverageExpression),Avg.MHCI.exprs=(re.avexprs$AverageExpression)$Avg.1)
pheno.re <- join(pheno.re, re.avg, type="left", match="all")
ti.avg <- data.frame(rnaseq.id=rownames(re.ti.avexprs$AverageExpression),Avg.MHCI.exprs=(re.ti.avexprs$AverageExpression)$Avg.1)
pheno.re.ti <- join(pheno.re.ti, ti.avg, type="left", match="all")
sc.avg <- data.frame(rnaseq.id=rownames(re.sc.avexprs$AverageExpression),Avg.MHCI.exprs=(re.sc.avexprs$AverageExpression)$Avg.1)
pheno.re.sc <- join(pheno.re.sc, sc.avg, type="left", match="all")
save(re,pheno.re,re.ti,pheno.re.ti,re.sc,pheno.re.sc,info,readme.re.norm,file="NORMALISED_COUNTS/RE_GRCh38_ByTissue_NormalisedCounts.RData")
write.table(pheno.re, "NORMALISED_COUNTS/RE_GRCh38_UndifferentiatedSamples_AverageMHCIExpression.txt", row.names=F, sep="\t", quote=F)

af.avg <- data.frame(rnaseq.id=rownames(af.avexprs$AverageExpression),Avg.MHCI.exprs=(af.avexprs$AverageExpression)$Avg.1)
pheno.af <- join(pheno.af, af.avg, type="left", match="all")
ti.avg <- data.frame(rnaseq.id=rownames(af.ti.avexprs$AverageExpression),Avg.MHCI.exprs=(af.ti.avexprs$AverageExpression)$Avg.1)
pheno.af.ti <- join(pheno.af.ti, ti.avg, type="left", match="all")
sc.avg <- data.frame(rnaseq.id=rownames(af.sc.avexprs$AverageExpression),Avg.MHCI.exprs=(af.sc.avexprs$AverageExpression)$Avg.1)
pheno.af.sc <- join(pheno.af.sc, sc.avg, type="left", match="all")
duo.avg <- data.frame(rnaseq.id=rownames(af.duo.avexprs$AverageExpression),Avg.MHCI.exprs=(af.duo.avexprs$AverageExpression)$Avg.1)
pheno.af.duo <- join(pheno.af.duo, duo.avg, type="left", match="all")
save(af,pheno.af,af.ti,pheno.af.ti,af.sc,pheno.af.sc,af.duo,pheno.af.duo,info,readme.af.norm,file="NORMALISED_COUNTS/AF_GRCh38_ByTissue_NormalisedCounts.RData")
write.table(pheno.af, "NORMALISED_COUNTS/AF_GRCh38_AverageMHCIExpression.txt", row.names=F, sep="\t", quote=F)


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


# AF DATASETS:

a.pca <- prcomp(t(af))
a.pca <- as.data.frame(a.pca$x)[,1:10]
a.pca$rnaseq.id <- rownames(a.pca)
a.pca <- join(pheno.af, a.pca, type="left", match="all")
write.table(a.pca, "NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_ALL_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_ALL_NormalisedCounts_PC1vPC2.pdf", width=15, height=5)
p1 <- ggplot(a.pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(a.pca, aes(PC1, PC2, fill=tissue))
p3 <- ggplot(a.pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite2,p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=3)
dev.off()

pca <- prcomp(t(af.ti))
pca <- as.data.frame(pca$x)[,1:10]
pca$rnaseq.id <- rownames(pca)
pca <- join(pheno.af.ti, pca, type="left", match="all")
write.table(pca, "NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_TI_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_TI_NormalisedCounts_PC1vPC2.pdf", width=15, height=5)
p1 <- ggplot(pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(pca, aes(PC1, PC2, fill=treatment))
p3 <- ggplot(pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_treat,p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=3)
dev.off()

pca <- prcomp(t(af.sc))
pca <- as.data.frame(pca$x)[,1:10]
pca$rnaseq.id <- rownames(pca)
pca <- join(pheno.af.sc, pca, type="left", match="all")
write.table(pca, "NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_SC_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_SC_NormalisedCounts_PC1vPC2.pdf", width=15, height=5)
p1 <- ggplot(pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(pca, aes(PC1, PC2, fill=treatment))
p3 <- ggplot(pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_treat,p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=3)
dev.off()

pca <- prcomp(t(af.duo))
pca <- as.data.frame(pca$x)[,1:10]
pca$rnaseq.id <- rownames(pca)
pca <- join(pheno.af.duo, pca, type="left", match="all")
write.table(pca, "NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_DUO_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_DUO_NormalisedCounts_PC1vPC2.pdf", width=10, height=10)
p1 <- ggplot(pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(pca, aes(PC1, PC2, fill=treatment))
p3 <- ggplot(pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_treat,p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=3)
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


# AF DATASETS:

pca <- af[(apply(af[,c(1:ncol(af))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.af[c("sex","diagnosis","treatment","tissue")]
meta_continuous <- pheno.af[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis","Treatment","Tissue")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_ALL_NormalisedCounts_HeatScree.pdf", width=10, height=7)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

pca <- af.ti[(apply(af.ti[,c(1:ncol(af.ti))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.af.ti[c("sex","diagnosis","treatment")]
meta_continuous <- pheno.af.ti[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis","Treatment")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_TI_NormalisedCounts_HeatScree.pdf", width=10, height=7)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

pca <- af.sc[(apply(af.sc[,c(1:ncol(af.sc))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.af.sc[c("sex","diagnosis","treatment")]
meta_continuous <- pheno.af.sc[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis","Treatment")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_SC_NormalisedCounts_HeatScree.pdf", width=10, height=7)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

pca <- af.duo[(apply(af.duo[,c(1:ncol(af.duo))],1,median))>=5,]
pca <- prcomp(t(pca))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- pheno.af.duo[c("sex","diagnosis","treatment")]
meta_continuous <- pheno.af.duo[c("age")]
colnames(meta_categorical) <- c("Sex","Diagnosis","Treatment")
colnames(meta_continuous) <- c("Age")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_DUO_NormalisedCounts_HeatScree.pdf", width=10, height=7)
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


# AF DATASETS:

norm.sex <- subset(af.ti, rownames(af.ti) %in% sex.genes)

Diagnosis <- c("skyblue","grey75")
Diagnosis <- Diagnosis[as.numeric(as.factor(pheno.af.ti$diagnosis))]
Sex <- c("darkred","darkblue")
Sex <- Sex[as.numeric(as.factor(pheno.af.ti$sex))]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

pdf("NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_TI_SexLinkedGenes_EuclideanAverageClustering.pdf", width=20, height=11)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=cbind(Diagnosis,Sex), dend=dend, rowLabels=c("Diagnosis","Sex"))
dev.off()

norm.sex <- subset(af.sc, rownames(af.sc) %in% sex.genes)

Diagnosis <- c("skyblue","grey75")
Diagnosis <- Diagnosis[as.numeric(as.factor(pheno.af.sc$diagnosis))]
Sex <- c("darkred","darkblue")
Sex <- Sex[as.numeric(as.factor(pheno.af.sc$sex))]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

pdf("NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_SC_SexLinkedGenes_EuclideanAverageClustering.pdf", width=20, height=11)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=cbind(Diagnosis,Sex), dend=dend, rowLabels=c("Diagnosis","Sex"))
dev.off()

norm.sex <- subset(af.duo, rownames(af.duo) %in% sex.genes)

Diagnosis <- c("skyblue","grey75")
Diagnosis <- Diagnosis[as.numeric(as.factor(pheno.af.duo$diagnosis))]
Sex <- c("darkred","darkblue")
Sex <- Sex[as.numeric(as.factor(pheno.af.duo$sex))]

d <- dist(t(norm.sex))
hc <- hclust(d, method="average")
dend <- as.dendrogram(hc)

pdf("NORMALISED_COUNTS_QC/AF_GRCh38/AF_GRCh38_DUO_SexLinkedGenes_EuclideanAverageClustering.pdf", width=20, height=11)
par(mar=c(10,7,2,4))
plot(dend)
colored_bars(colors=cbind(Diagnosis,Sex), dend=dend, rowLabels=c("Diagnosis","Sex"))
dev.off()


################
# SESSION INFO #
################

sessionInfo()
rm(list=ls())

