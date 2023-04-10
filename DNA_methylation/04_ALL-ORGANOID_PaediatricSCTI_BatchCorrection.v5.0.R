# Batch correct all organoid data using ComBat, including the "METH-K.NAYAK-20038_17JUN22" batch
# Updated to subset the "all organoid" dataset to just paediatric SC & TI samples, dropping any treated samples except those belonging to the AZA experiment
# Updated from 03_EPIC_BatchCorrection.v4.2.R to be more in line with Rachel Edgar's NLRC5/MHCI methylation array analysis scripts
# 11Jul22
#
# Folders required:
#     METHYLATION_ARRAY_ANALYSIS/SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA
#
# Files required:
#     METHYLATION_ARRAY_ANALYSIS/R_DATASETS/AllOrganoids_FilteredNormalisedNoBatchCorrection_beta.RData (containing normalised, QCed dataset for all combined K450 & EPIC organoid samples produced by 02_ALL-ORGANOID_K450-EPIC_CombineData.v3.0.R)
#
# Useful user guides/tutorials:
#     sva: https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf


#######################
# WORKING ENVIRONMENT #
#######################

setwd("/home/fp215/rds/rds-fp215-working/METHYLATION_ARRAY_ANALYSIS")

library(genefilter)
library(matrixStats)
library(minfi)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(IlluminaHumanMethylationEPICmanifest)
library(sva)
library(plyr)
library(lattice)
library(ggplot2)
library(rafalib)
library(reshape)
library(ClassDiscovery)
source("/home/fp215/SCRIPTS/ADDENBROOKES_SCRIPTS/GENERAL_PURPOSE/PlotTools.v3.1.R")

PCs_to_view <- 10

# ALL ORGANOID SAMPLES:

load("R_DATASETS/AllOrganoids_FilteredNormalisedNoBatchCorrection_beta.RData")

readme.beta

samples$sentrix.id <- as.character(samples$sentrix.id)

# SUBSET TO JUST PAEDIATRIC SC & TI SAMPLES:
# Also remove any treated samples that aren't part of the new AZA experiment & any "NK" or "Other.GI" samples

select <- c("NT","AZA-3W","AZA-72H")

samples <- subset(samples, age.group=="paediatric" & (tissue=="TI" | tissue=="SC"))
samples <- subset(samples, treatment.1 %in% select)
samples <- subset(samples, diagnosis!="NK" & diagnosis!="Other.GI")

overlap.beta <- overlap.beta[,samples$array.id]

dim(overlap.beta)
nrow(samples)

# Pull out IBD & control only datasets:

ibd.samples <- subset(samples, diagnosis=="CD" | diagnosis=="UC" | diagnosis=="IBD-U")
cont.samples <- subset(samples, diagnosis=="Control")

# Load pre-saved colour scheme for the sample batches:

load("R_DATASETS/ArrayColourSchemes.RData")

readme.batches

all.batches <- subset(all.batches, Batch %in% samples$sentrix.id)
ibd.batches <- subset(all.batches, Batch %in% ibd.samples$sentrix.id)
cont.batches <- subset(all.batches, Batch %in% cont.samples$sentrix.id)


#############
# FUNCTIONS #
#############

# Impute NAs (with row mean), 0s & 1s so that can use M values for batch correction:

imputeMedian <- function(x) apply(x, 1, function(x){x[is.na(x)] <- median(x, na.rm=T); x})

# Convert beta values to M values:

convert <- function(x){log2(x/(1-x))}

# Convert M values to beta values:

reconvert <- function(x){2^x/((2^x)+1)}


######################################################
# ALL ORGANOID SAMPLES: BATCH CORRECTION, RE VERSION #
######################################################

# EDIT THE BETA VALUES FOR COMBAT #1:

orgi.beta <- overlap.beta
orgi.beta[orgi.beta==0] <- 0.01
orgi.beta[orgi.beta==1] <- 0.99
orgi.beta <- t(imputeMedian(orgi.beta))

# EXTRACT M VALUES #1:

orgi.mval <- structure(sapply(orgi.beta, convert), dim=dim(orgi.beta))
rownames(orgi.mval) <- rownames(orgi.beta)
colnames(orgi.mval) <- colnames(orgi.beta)

# BATCH CORRECTION #1:

orgi.mval.re1 <- ComBat(dat=orgi.mval, batch=samples$array.type, mod=NULL, par.prior=T)

# CONVERT BACK TO BETA VALUES #1:

orgi.beta.re1 <- structure(sapply(orgi.mval.re1, reconvert), dim=dim(orgi.mval.re1))
rownames(orgi.beta.re1) <- rownames(orgi.mval.re1)
colnames(orgi.beta.re1) <- colnames(orgi.mval.re1)

# HEAT SCREE AFTER COMBAT FOR ARRAY TYPE #1:

pca <- prcomp(t(orgi.beta.re1))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- samples[c("diagnosis","sex","tissue","passage.or.rescope.no","sentrix.id","RE.batch","array.type")]
meta_continuous <- samples[c("age")]
colnames(meta_categorical) <- c("Diagnosis","Sex","Tissue","Passage","Sentrix ID","Batch","Array")
colnames(meta_continuous) <- "Age"
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA/AllOrganoids-PaediatricSCTI_BatchCorrection-1-ArrayType_HeatScreePlot.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

# SEPARATE OUT THE IBD & CONTROL SAMPLES #1:

ibd.beta <- orgi.beta.re1[,ibd.samples$array.id]
cont.beta <- orgi.beta.re1[,cont.samples$array.id]

# EDIT THE BETA VALUES FOR COMBAT #2:

orgi.beta <- orgi.beta.re1
orgi.beta[orgi.beta==0] <- 0.01
orgi.beta[orgi.beta==1] <- 0.99
orgi.beta <- t(imputeMedian(orgi.beta))

# EXTRACT M VALUES #2:

orgi.mval <- structure(sapply(orgi.beta, convert), dim=dim(orgi.beta))
rownames(orgi.mval) <- rownames(orgi.beta)
colnames(orgi.mval) <- colnames(orgi.beta)

# BATCH CORRECTION #2:

mod <- model.matrix(~as.factor(sample.id), data=samples)
orgi.mval.re2 <- ComBat(dat=orgi.mval, batch=samples$RE.batch, mod=mod, par.prior=T)

# CONVERT BACK TO BETA VALUES #2:

orgi.beta.re2 <- structure(sapply(orgi.mval.re2, reconvert), dim=dim(orgi.mval.re2))
rownames(orgi.beta.re2) <- rownames(orgi.mval.re2)
colnames(orgi.beta.re2) <- colnames(orgi.mval.re2)

# HEAT SCREE AFTER COMBAT FOR ARRAY TYPE #2:

pca <- prcomp(t(orgi.beta.re2))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- samples[c("diagnosis","sex","tissue","passage.or.rescope.no","sentrix.id","RE.batch","array.type")]
meta_continuous <- samples[c("age")]
colnames(meta_categorical) <- c("Diagnosis","Sex","Tissue","Passage","Sentrix ID","Batch","Array")
colnames(meta_continuous) <- "Age"
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA/AllOrganoids-PaediatricSCTI_BatchCorrection-2-REBatch_HeatScreePlot.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

# SEPARATE OUT THE IBD & CONTROL SAMPLES #1:

ibd.beta <- orgi.beta.re2[,ibd.samples$array.id]
cont.beta <- orgi.beta.re2[,cont.samples$array.id]

# PCA PLOTS AFTER COMBAT FOR ARRAY TYPE #2:
# NB, no K450 array IBD samples

all.pca <- prcomp(t(orgi.beta.re2))
all.pca <- as.data.frame(all.pca$x)
all.pca$array.id <- rownames(all.pca)
all.pca <- join(all.pca, samples, type="left", match="all")
write.table(all.pca, "SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA/AllOrganoids-PaediatricSCTI_ALL_BatchCorrection-2-REBatch_PCs.txt", row.names=F, sep="\t", quote=F)

ibd.pca <- prcomp(t(ibd.beta))
ibd.pca <- as.data.frame(ibd.pca$x)
ibd.pca$array.id <- rownames(ibd.pca)
ibd.pca <- join(ibd.pca, ibd.samples, type="left", match="all")
write.table(ibd.pca, "SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA/AllOrganoids-PaediatricSCTI_IBD_BatchCorrection-2-REBatch_PCs.txt", row.names=F, sep="\t", quote=F)

cont.pca <- prcomp(t(cont.beta))
cont.pca <- as.data.frame(cont.pca$x)
cont.pca$array.id <- rownames(cont.pca)
cont.pca <- join(cont.pca, samples, type="left", match="all")
write.table(cont.pca, "SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA/AllOrganoids-PaediatricSCTI_Controls_BatchCorrection-2-REBatch_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA/AllOrganoids-PaediatricSCTI_BatchCorrection-2-REBatch_PC1vPC2.pdf", width=15, height=16)
p1 <- ggplot(all.pca, aes(PC1, PC2, fill=diagnosis, color=tissue))
p2 <- ggplot(ibd.pca, aes(PC1, PC2, fill=diagnosis, color=tissue))
p3 <- ggplot(cont.pca, aes(PC1, PC2, fill=diagnosis, color=tissue))
p4 <- ggplot(all.pca, aes(PC1, PC2, fill=tissue))
p5 <- ggplot(ibd.pca, aes(PC1, PC2, fill=tissue))
p6 <- ggplot(cont.pca, aes(PC1, PC2, fill=tissue))
p7 <- ggplot(all.pca, aes(PC1, PC2))
p8 <- ggplot(ibd.pca, aes(PC1, PC2))
p9 <- ggplot(cont.pca, aes(PC1, PC2))
grid.arrange(p1+geom_point(shape=21,size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2+colscale_sampsite3+ggtitle("ALL")+theme(plot.title=element_text(face="bold")),
             p2+geom_point(shape=21,size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2+colscale_sampsite3+ggtitle("IBD")+theme(plot.title=element_text(face="bold")),
             p3+geom_point(shape=21,size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2+colscale_sampsite3+ggtitle("Control")+theme(plot.title=element_text(face="bold")),
             p4+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite,
             p5+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite,
             p6+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite,
             p7+geom_point(aes(shape=RE.batch, fill=array.type, color="black"),size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_array+colscale_array+scale_shape_manual(name="Batch",values=c(22,24,21),labels=c("NEW","OLD","RE")),
             p8+geom_point(aes(shape=RE.batch, fill=array.type, color="black"),size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_array+colscale_array+scale_shape_manual(name="Batch",values=c(22,24,21),labels=c("NEW","OLD","RE")),
             p9+geom_point(aes(shape=RE.batch, fill=array.type, color="black"),size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_array+colscale_array+scale_shape_manual(name="Batch",values=c(22,24,21),labels=c("NEW","OLD","RE")),ncol=3)
dev.off()

# SAVE BATCH CORRECTED DATASETS:

readme.combat.re <- data.frame(Dataset=c("overlap.beta","orgi.beta.re1","orgi.mval.re1","orgi.beta.re2","orgi.mval.re2","samples"),Description=c("Funnorm normalised, filtered, K450+EPIC paediatric organoid SC & TI beta values","overlap.beta ComBat corrected for array type (beta values)","overlap.beta ComBat corrected for array type (M values)","overlap.beta.re1 ComBat corrected for RE.batch (beta values)","overlap.beta.re1 ComBat corrected for RE.batch (M values)","Sample Information"))
save(overlap.beta,orgi.beta.re1,orgi.mval.re1,orgi.beta.re2,orgi.mval.re2,samples,readme.combat.re,file="R_DATASETS/AllOrganoids-PaediatricSCTI_RE-BatchCorrected.RData")


###################################################
# ALL ORGANOID SAMPLES: BATCH CORRECTION, BY CHIP #
###################################################

# BATCH CORRECTION:

mod <- model.matrix(~as.factor(tissue), data=samples)
orgi.beta.chip <- ComBat(dat=overlap.beta, batch=samples$sentrix.id, mod=mod)

# EXTRACT M VALUES:

orgi.mval.chip <- structure(sapply(orgi.beta.chip, convert), dim=dim(orgi.beta.chip))
rownames(orgi.mval.chip) <- rownames(orgi.beta.chip)
colnames(orgi.mval.chip) <- colnames(orgi.beta.chip)

# HEAT SCREE AFTER COMBAT FOR ARRAY TYPE:

pca <- prcomp(t(orgi.beta.chip))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- samples[c("diagnosis","sex","tissue","passage.or.rescope.no","sentrix.id","RE.batch","array.type")]
meta_continuous <- samples[c("age")]
colnames(meta_categorical) <- c("Diagnosis","Sex","Tissue","Passage","Sentrix ID","Batch","Array")
colnames(meta_continuous) <- "Age"
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA/AllOrganoids-PaediatricSCTI_BatchCorrection-SentrixID_HeatScreePlot.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

# SEPARATE OUT THE IBD & CONTROL SAMPLES:

ibd.beta <- orgi.beta.chip[,ibd.samples$array.id]
cont.beta <- orgi.beta.chip[,cont.samples$array.id]

# PCA PLOTS AFTER COMBAT FOR ARRAY TYPE:

all.pca <- prcomp(t(orgi.beta.chip))
all.pca <- as.data.frame(all.pca$x)
all.pca$array.id <- rownames(all.pca)
all.pca <- join(all.pca, samples, type="left", match="all")
write.table(all.pca, "SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA/AllOrganoids-PaediatricSCTI_ALL_BatchCorrection-SentrixID_PCs.txt", row.names=F, sep="\t", quote=F)

ibd.pca <- prcomp(t(ibd.beta))
ibd.pca <- as.data.frame(ibd.pca$x)
ibd.pca$array.id <- rownames(ibd.pca)
ibd.pca <- join(ibd.pca, ibd.samples, type="left", match="all")
write.table(ibd.pca, "SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA/AllOrganoids-PaediatricSCTI_IBD_BatchCorrection-SentrixID_PCs.txt", row.names=F, sep="\t", quote=F)

cont.pca <- prcomp(t(cont.beta))
cont.pca <- as.data.frame(cont.pca$x)
cont.pca$array.id <- rownames(cont.pca)
cont.pca <- join(cont.pca, samples, type="left", match="all")
write.table(cont.pca, "SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA/AllOrganoids-PaediatricSCTI_Controls_BatchCorrection-SentrixID_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("SAMPLE_QC_ALL-ORGANOID/07_BATCH_CORRECTED_DATA/AllOrganoids-PaediatricSCTI_BatchCorrection-SentrixID_PC1vPC2.pdf", width=15, height=21)
p1 <- ggplot(all.pca, aes(PC1, PC2, fill=diagnosis, color=tissue))
p2 <- ggplot(ibd.pca, aes(PC1, PC2, fill=diagnosis, color=tissue))
p3 <- ggplot(cont.pca, aes(PC1, PC2, fill=diagnosis, color=tissue))
p4 <- ggplot(all.pca, aes(PC1, PC2, fill=tissue))
p5 <- ggplot(ibd.pca, aes(PC1, PC2, fill=tissue))
p6 <- ggplot(cont.pca, aes(PC1, PC2, fill=tissue))
p7 <- ggplot(all.pca, aes(PC1, PC2))
p8 <- ggplot(ibd.pca, aes(PC1, PC2))
p9 <- ggplot(cont.pca, aes(PC1, PC2))
grid.arrange(p1+geom_point(shape=21,size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2+colscale_sampsite3+ggtitle("ALL")+theme(plot.title=element_text(face="bold")),
             p2+geom_point(shape=21,size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2+colscale_sampsite3+ggtitle("IBD")+theme(plot.title=element_text(face="bold")),
             p3+geom_point(shape=21,size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2+colscale_sampsite3+ggtitle("Control")+theme(plot.title=element_text(face="bold")),
             p4+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite,
             p5+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite,
             p6+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite,
             p7+geom_point(aes(shape=RE.batch, fill=array.type, color="black"),size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_array+colscale_array+scale_shape_manual(name="Batch",values=c(22,24,21),labels=c("NEW","OLD","RE")),
             p8+geom_point(aes(shape=RE.batch, fill=array.type, color="black"),size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_array+colscale_array+scale_shape_manual(name="Batch",values=c(22,24,21),labels=c("NEW","OLD","RE")),
             p9+geom_point(aes(shape=RE.batch, fill=array.type, color="black"),size=2)+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_array+colscale_array+scale_shape_manual(name="Batch",values=c(22,24,21),labels=c("NEW","OLD","RE")),ncol=3)
dev.off()

# SAVE BATCH CORRECTED DATASET:

readme.combat.chip <- data.frame(Dataset=c("overlap.beta","orgi.beta.chip","orgi.mval.chip","samples"),Description=c("Funnorm normalised, filtered, K450+EPIC paediatric organoid SC & TI beta values","overlap.beta ComBat corrected for Sentrix ID (beta values)","overlap.beta ComBat corrected for Sentrix ID (M values)","Sample Information"))
save(overlap.beta,orgi.beta.chip,orgi.mval.chip,samples,readme.combat.chip,file="R_DATASETS/AllOrganoids-PaediatricSCTI_SentrixID-BatchCorrected.RData")


################
# SESSION INFO #
################

sessionInfo()
rm(list=ls())