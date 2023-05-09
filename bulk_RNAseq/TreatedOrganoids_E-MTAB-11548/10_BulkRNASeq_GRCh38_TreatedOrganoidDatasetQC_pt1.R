# QC OF THE NGS-K.NAYAK-40706_RE (E-MTAB-11548) RAW BULK RNASEQ COUNTS DATASET - Part 1
# NOISeq & DESeq2 QC plots, check for batch effects
# Data aligned to GRCh38
# 24Octp22 (fp215)
# 
# Folders required:
#     - RAW_COUNTS_QC
#     - RAW_COUNTS_QC/RE_GRCh38
#     - NORMALISED_COUNTS_QC
#     - NORMALISED_COUNTS_QC/RE_GRCh38
#       
# Files required:
#     - RE raw counts: /home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/RE_GRCh38_per-gene-counts.2.txt
#     - Sample info: /home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/All_bulkRNASeq_UpdatedSamples_24Oct22.RData
#
# Useful documentation:
#     - biomaRt: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
#     - NOISeq documentation: https://bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
#     - DESeq2 documentation: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#     - SVA documentation: https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf
#
# Functions used:
#     - SV correlation function: /home/fp215/SCRIPTS/ADDENBROOKES_SCRIPTS/GENERAL_PURPOSE/BatchEffectCorrelation_function.v1.R
#     - Plot colour schemes & themes: /home/fp215/SCRIPTS/ADDENBROOKES_SCRIPTS/GENERAL_PURPOSE/PlotTools.v3.4.R


#######################
# WORKING ENVIRONMENT #
#######################

setwd("/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER")
library(biomaRt)
library(NOISeq)
library(DESeq2)
library(ggplot2)
library(ggpubr)
require(gridExtra)
library(pheatmap)
library(RColorBrewer)
library(sva)
library(plyr)
source("/home/fp215/SCRIPTS/ADDENBROOKES_SCRIPTS/GENERAL_PURPOSE/BatchEffectCorrelation_function.v1.R")
source("/home/fp215/SCRIPTS/ADDENBROOKES_SCRIPTS/GENERAL_PURPOSE/PlotTools.v3.4.R")

cols <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)


#####################################
# EXTRACT ENSEMBL GRCh38 ANNOTATION #
#####################################

# Extract annotation for GRCh38 release 99:

ensembl <- useEnsembl(biomart="ensembl", version=99, dataset="hsapiens_gene_ensembl")

# Pull out information of interest:

chrs <- c(1:22,"X","Y")
info <- getBM(attributes=c("external_gene_name", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype", "transcript_count"), mart=ensembl)
info$length <- info$end_position - info$start_position
info <- rename(info, c("external_gene_name"="gene.symbol","ensembl_gene_id"="gene","chromosome_name"="chr","start_position"="start","end_position"="end","percentage_gene_gc_content"="percent.gc","gene_biotype"="biotype","transcript_count"="tc"))
info <- subset(info, chr %in% chrs)

# Save for future use:

save(info,file="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/GRCh38_release-99_AnnotationInfo.RData")


#############
# LOAD DATA #
#############

# Raw counts (including RE differentiated samples):

raw.re.all <- read.table("/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/RE_GRCh38_per-gene-counts.2.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(raw.re.all) <- raw.re.all$Geneid
raw.re.all$Geneid <- NULL

# Sample info:

load("All_bulkRNASeq_UpdatedSamples_24Oct22.RData")

readme.bulk

samples <- join(bulk.sample, bulk.clin, type="left", match="all")
samples$differentiated[is.na(samples$differentiated)] <- "undifferentiated"

table(samples$differentiated)

pheno.re.all <- data.frame(rnaseq.id=colnames(raw.re.all))
pheno.re.all <- join(pheno.re.all, samples, type="left", match="all")


#################################
# REMOVE DIFFERENTIATED SAMPLES #
#################################

pheno.re <- subset(pheno.re.all, differentiated=="undifferentiated")
raw.re <- raw.re.all[,pheno.re$rnaseq.id]


#####################
# SAVE RAW DATASETS #
#####################

readme.re.raw <- data.frame(Dataset=c("raw.re.all","pheno.re.all","raw.re","pheno.re","info"),Description=c("NGS-K.Nayak-40706 (RE) raw counts including differentiated samples (n=42)","Corresponding sample info, updated 24Oct22","NGS-K.Nayak-40706 (RE) raw counts excluding differentiated samples (n=30)","Corresponding samlpe linfo, updated 24Oct22","Ensembl GRCh38 release 99 gene annotation info"))
save(raw.re.all,pheno.re.all,raw.re,pheno.re,info,readme.re.raw,file="RAW_COUNTS/RE_GRCh38_RawCounts.RData")


#####################################
# EXTRACT NOISEQ QC ANNOTATION INFO #
#####################################

length <- unique(info[c("gene","length")])
rownames(length) <- length[,1]
length$gene <- NULL

gc <- unique(info[c("gene","percent.gc")])
rownames(gc) <- gc[,1]
gc$gene <- NULL

biotype <- unique(info[c("gene","biotype")])
rownames(biotype) <- biotype[,1]
biotype$gene <- NULL

coords <- unique(info[c("gene","chr","start","end")])
rownames(coords) <- coords[,1]
coords$gene <- NULL


############################
# CREATE NOISEQ QC DATASET #
############################

re.noiseq <- readData(data=raw.re, length=length, gc=gc, biotype=biotype, chromosome=coords, factors=pheno.re)


###################
# NOISEQ QC PLOTS #
###################

# BIOTYPE DISTRIBUTION:
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/Biodetection

biodetection <- dat(re.noiseq, k=0, type="biodetection", factor="tissue")

pdf("RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_BiotypeDistribution.pdf", width=10, height=5)
explo.plot(biodetection, samples="TI", plottype="persample", toplot="protein_coding")
explo.plot(biodetection, samples="SC", plottype="persample", toplot="protein_coding")
dev.off()

# COUNT DISTRIBUTION:
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/CountsBio

counts <- dat(re.noiseq, factor="tissue", type="countsbio")

pdf("RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_CountDistribution.pdf", width=10, height=6)
explo.plot(counts, samples="TI", toplot="global", plottype="boxplot")
explo.plot(counts, samples="SC", toplot="global", plottype="boxplot")
dev.off()

# PROTEIN CODING COUNT DISTRIBUTION:

r.counts <- dat(re.noiseq, factor=NULL, type="countsbio")

pdf("RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_ProteinCoding_CountDistribution.pdf", width=30, height=5)
explo.plot(r.counts, toplot="protein_coding", samples=NULL, plottype="boxplot")
dev.off()

# PERCENTAGE LOW COUNT FEATURES PER SAMPLE:

pdf("RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_ProteinCoding_PercentLowCount.pdf", width=30, height=5)
explo.plot(r.counts, toplot="protein_coding", samples=NULL, plottype="barplot")
dev.off()

# CHECK FOR LENGTH BIAS (significant p-value + R2 > 70%):
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/lengthbias

len <- dat(re.noiseq, factor="tissue", type="lengthbias")

pdf("RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_ProteinCoding_LengthBias.pdf", width=7, height=5)
explo.plot(len, samples=NULL, toplot="protein_coding")
dev.off()

# CHECK FOR GC BIAS (significant p-value + R2 > 70%):
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/GCbias

raw.gc <- dat(re.noiseq, factor="tissue", type="GCbias")

pdf("RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_ProteinCoding_GCBias.pdf", width=7, height=5)
explo.plot(raw.gc, samples=NULL, toplot="protein_coding")
dev.off()


#############################
# CREATE DESEQ2 QC DATASETS #
#############################

re.factors <- pheno.re
re.factors$group <- paste0(re.factors$tissue,".",re.factors$diagnosis)
re.factors$group2 <- paste0(re.factors$group,".",re.factors$treatment)
re.factors$group3 <- paste0(re.factors$diagnosis,".",re.factors$treatment)
rownames(re.factors) <- re.factors$rnaseq.id
re.factors$diagnosis <- as.factor(re.factors$diagnosis)

af.dds <- DESeqDataSetFromMatrix(countData=raw.af, colData=af.factors, design=~group2)
af.dds$diagnosis <- relevel(af.dds$diagnosis, ref="Control")
af.dds <- af.dds[rowSums(counts(af.dds))>1,]
af.vst <- vst(af.dds,blind=F)


###################
# DESEQ2 QC PLOTS #
###################

# EUCLDIDEAN DISTANCE HEATMAP:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-sample-to-sample-distances

raw.dist <- dist(t(assay(re.vst)))

pdf("RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_RawEuclideanDistanceHeatmap.pdf", width=20, height=20)
matrix <- as.matrix(raw.dist)
rownames(matrix) <- re.vst$group3
colnames(matrix) <- re.vst$tissue
pheatmap(matrix, clustering_distance_rows=raw.dist, clustering_distance_cols=raw.dist, col=cols, cellwidth=10, cellheight=10)
dev.off()

# COOK'S CUTOFF FOR OUTLIERS:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#approach-to-count-outliers

r.dds <- DESeq(re.dds)
save(re.dds,re.vst,r.dds,pheno.re,file="RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_RawDataQCInput.RData")

pdf("RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_CooksCutoff_RawNoBatchCorrection.pdf", width=30, height=7)
boxplot(log10(assays(r.dds)[["cooks"]]), range=0, las=2)
dev.off()

# PCA PLOTS:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#principal-component-plot-of-the-samples

re <- counts(r.dds, normalized=T)

readme.re.norm <- data.frame(Dataset=c("re","pheno.re","info"),Description=c("NGS-K.Nayak-40706 (RE) normalised counts (n=30)","RE sample info, updated 24Oct22","Ensembl GRCh37 release 107 gene annotation info"))
save(re,pheno.re,info,readme.re.norm,file="NORMALISED_COUNTS/RE_GRCh38_NormalisedCounts.RData")

pca.raw <- plotPCA(re.vst, intgroup=c("group2"), returnData=T)
percentvar <- round(100*attr(pca.raw, "percentVar"))

pdf("RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_Raw_Tissue-Diagnosis-Treatment_PC1vPC2.pdf", width=10, height=10)
ggplot(pca.raw, aes(PC1, PC2, color=group2))+geom_point(size=3)+xlab(paste0("PC1: ",percentvar[1],"% variance"))+ylab(paste0("PC2: ",percentvar[2],"%variance"))
dev.off()

r.pca <- prcomp(t(re))
r.pca <- as.data.frame(r.pca$x)[,1:10]
r.pca$rnaseq.id <- rownames(r.pca)
r.pca <- join(pheno.re, r.pca, type="left", match="all")
write.table(r.pca, "NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_ALL_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/RE_GRCh38/RE_GRCh38_ALL_NormalisedCounts_PC1vPC2.pdf", width=10, height=10)
p1 <- ggplot(r.pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(r.pca, aes(PC1, PC2, fill=tissue))
p3 <- ggplot(r.pca, aes(PC1, PC2, fill=treatment))
p4 <- ggplot(r.pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite2,p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_treat,p4+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=2)
dev.off()


###########################
# CHECK FOR BATCH EFFECTS #
###########################

# Create a matrix of the normalised counts and a dataframe of the factors for sva:

r.norm <- counts(r.dds, normalized=T)
r.svafactors <- as.data.frame(re.factors)

# Create models (mod=full model, mod0=null, NB, including known batch effects: diagnosis, tissue & treatment):

r.mod <- model.matrix(~as.factor(group2), data=r.svafactors)
r.mod0 <- model.matrix(~1, data=r.svafactors)

# Calculate the first 5 SVs:

r.svseq <- svaseq(r.norm,r.mod,r.mod0,n.sv=5)

# Plot variables against the predicted surrogate variables:

meta_categorical <- data.frame(tissue=pheno.re$tissue,diagnosis=pheno.re$diagnosis,sex=pheno.re$sex,treatment=pheno.re$treatment)
meta_continuous <- data.frame(age=pheno.re$age,passage=as.numeric(gsub("P","",pheno.re$passage)))
meta_sva <- data.frame(SV1=r.svseq$sv[,1],SV2=r.svseq$sv[,2],SV3=r.svseq$sv[,3],SV4=r.svseq$sv[,4],SV5=r.svseq$sv[,5])

r.cor <- SVcor(meta_categorical,meta_continuous,meta_sva)

pdf("RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_SVACorrelations_5SVs.pdf", width=8, height=7)
r.cor$plot.input
dev.off()

write.table(r.cor$sv.correlation, "RAW_COUNTS_QC/RE_GRCh38/RE_GRCh38_SVACorrelations_5SVs.txt", row.names=F, sep="\t", quote=F)


################
# SESSION INFO #
################

sessionInfo()
rm(list=ls())

