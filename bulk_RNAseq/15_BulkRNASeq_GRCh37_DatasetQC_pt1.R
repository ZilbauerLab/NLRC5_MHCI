# QC OF THE E-MTAB-5464 (HOWELL) & NGS-K.NAYAK-40706_RE + NGS-K.NAYAK-40750_AF (RE-AF) RAW BULK RNASEQ COUNTS DATASETS - Part 1
# NOISeq & DESeq2 QC plots, check for batch effects
# Data aligned to GRCh37
# 30Sep22 (fp215)
# 24Oct22 (fp215) - updates using separated AF & RE datasets & updated sample info, also remove the differentiated samples in the RE dataset from the start
# 
# Folders required:
#     - RAW_COUNTS_QC
#     - RAW_COUNTS_QC/HOWELL_GRCh37
#     - RAW_COUNTS_QC/RE_GRCh37
#     - RAW_COUNTS_QC/AF_GRCh37
#     - NORMALISED_COUNTS_QC
#     - NORMALISED_COUNTS_QC/HOWELL_GRCh37
#     - NORMALISED_COUNTS_QC/RE_GRCh37
#     - NORMALISED_COUNTS_QC/AF_GRCh37
#       
# Files required:
#     - Howell raw counts: /home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/Howell_GRCh37_per-gene-counts.2.txt
#     - RE raw counts: /home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/RE_GRCh37_per-gene-counts.2.txt
#     - AF raw counts: /home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/AF_GRCh37_per-gene-counts.2.txt
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
# EXTRACT ENSEMBL GRCh37 ANNOTATION #
#####################################

# Extract annotation for GRCh37 release 107:

listEnsembl(GRCh=37)

ensembl <- useEnsembl(biomart="ensembl", GRCh=37, dataset="hsapiens_gene_ensembl")

# Pull out information of interest:

chrs <- c(1:22,"X","Y")
info <- getBM(attributes=c("external_gene_name", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype", "transcript_count"), mart=ensembl)
info$length <- info$end_position - info$start_position
info <- rename(info, c("external_gene_name"="gene.symbol","ensembl_gene_id"="gene","chromosome_name"="chr","start_position"="start","end_position"="end","percentage_gene_gc_content"="percent.gc","gene_biotype"="biotype","transcript_count"="tc"))
info <- subset(info, chr %in% chrs)

# Save for future use:

save(info,file="/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/GRCh37_release-107_AnnotationInfo.RData")


#############
# LOAD DATA #
#############

# Raw counts (including RE differentiated samples):

#raw.howell <- read.table("/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/Howell_GRCh37_per-gene-counts.2.txt", header=T, sep="\t", stringsAsFactors=F)
#rownames(raw.howell) <- raw.howell$Geneid
#raw.howell$Geneid <- NULL
raw.re.all <- read.table("/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/RE_GRCh37_per-gene-counts.2.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(raw.re.all) <- raw.re.all$Geneid
raw.re.all$Geneid <- NULL
raw.af <- read.table("/home/fp215/rds/rds-fp215-working/RNASEQ/NLRC5_PAPER/RAW_COUNTS/AF_GRCh37_per-gene-counts.2.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(raw.af) <- raw.af$Geneid
raw.af$Geneid <- NULL

# Sample info:

load("All_bulkRNASeq_UpdatedSamples_24Oct22.RData")

readme.bulk

samples <- join(bulk.sample, bulk.clin, type="left", match="all")
samples$differentiated[is.na(samples$differentiated)] <- "undifferentiated"

table(samples$differentiated)

#pheno.howell <- data.frame(rnaseq.id=colnames(raw.howell))
#pheno.howell <- join(pheno.howell, samples, type="left", match="all")

pheno.re.all <- data.frame(rnaseq.id=colnames(raw.re.all))
pheno.re.all <- join(pheno.re.all, samples, type="left", match="all")

pheno.af <- data.frame(rnaseq.id=colnames(raw.af))
pheno.af <- join(pheno.af, samples, type="left", match="all")


#################################
# REMOVE DIFFERENTIATED SAMPLES #
#################################

pheno.re <- subset(pheno.re.all, differentiated=="undifferentiated")
raw.re <- raw.re.all[,pheno.re$rnaseq.id]


#####################
# SAVE RAW DATASETS #
#####################

#readme.howell.raw <- data.frame(Dataset=c("raw.howell","pheno.howell","info"),Description=c("E-MTAB-5464 (Howell) raw counts (n=76)","Corresponding sample info, updated 29Sep22","Ensembl GRCh37 release 107 gene annotation info"))
#save(raw.howell,pheno.howell,info,readme.howell.raw,file="RAW_COUNTS/Howell_GRCh37_RawCounts.RData")

readme.re.raw <- data.frame(Dataset=c("raw.re.all","pheno.re.all","raw.re","pheno.re","info"),Description=c("NGS-K.Nayak-40706 (RE) raw counts including differentiated samples (n=42)","Corresponding sample info, updated 24Oct22","NGS-K.Nayak-40706 (RE) raw counts excluding differentiated samples (n=30)","Corresponding samlpe linfo, updated 24Oct22","Ensembl GRCh37 release 107 gene annotation info"))
save(raw.re.all,pheno.re.all,raw.re,pheno.re,info,readme.re.raw,file="RAW_COUNTS/RE_GRCh37_RawCounts.RData")
readme.af.raw <- data.frame(Dataset=c("raw.af","pheno.af","info"),Description=c("NGS-K.Nayak-40750 (AF) raw counts (n=90)","Corresponding sample info, updated 24Oct22","Ensembl GRCh37 release 107 gene annotation info"))
save(raw.af,pheno.af,info,readme.af.raw,file="RAW_COUNTS/AF_GRCh37_RawCounts.RData")


#####################################
# EXTRACT NOISEQ QC ANNOTATION INFO #
#####################################

enslength <- unique(info[c("gene","length")])
rownames(enslength) <- enslength[,1]
enslength$gene <- NULL
head(enslength,n=1)

ensgc <- unique(info[c("gene","percent.gc")])
rownames(ensgc) <- ensgc[,1]
ensgc$gene <- NULL
head(ensgc,n=1)

ensbiotype <- unique(info[c("gene","biotype")])
rownames(ensbiotype) <- ensbiotype[,1]
ensbiotype$gene <- NULL
head(ensbiotype,n=1)

enscoords <- unique(info[c("gene","chr","start","end")])
rownames(enscoords) <- enscoords[,1]
enscoords$gene <- NULL
head(enscoords,n=1)


#############################
# CREATE NOISEQ QC DATASETS #
#############################

#howell.noiseq <- readData(data=raw.howell, length=enslength, gc=ensgc, biotype=ensbiotype, chromosome=enscoords, factors=pheno.howell)
re.noiseq <- readData(data=raw.re, length=enslength, gc=ensgc, biotype=ensbiotype, chromosome=enscoords, factors=pheno.re)
af.noiseq <- readData(data=raw.af, length=enslength, gc=ensgc, biotype=ensbiotype, chromosome=enscoords, factors=pheno.af)


###################
# NOISEQ QC PLOTS #
###################

# BIOTYPE DISTRIBUTION:
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/Biodetection

#biodetection <- dat(howell.noiseq, k=0, type="biodetection", factor="tissue")

#pdf("RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_BiotypeDistribution.pdf", width=10, height=5)
#explo.plot(biodetection, samples="TI", plottype="persample", toplot="protein_coding")
#explo.plot(biodetection, samples="SC", plottype="persample", toplot="protein_coding")
#explo.plot(biodetection, samples="AC", plottype="persample", toplot="protein_coding")
#dev.off()

biodetection <- dat(re.noiseq, k=0, type="biodetection", factor="tissue")

pdf("RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_BiotypeDistribution.pdf", width=10, height=5)
explo.plot(biodetection, samples="TI", plottype="persample", toplot="protein_coding")
explo.plot(biodetection, samples="SC", plottype="persample", toplot="protein_coding")
dev.off()

biodetection <- dat(af.noiseq, k=0, type="biodetection", factor="tissue")

pdf("RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_BiotypeDistribution.pdf", width=10, height=5)
explo.plot(biodetection, samples="TI", plottype="persample", toplot="protein_coding")
explo.plot(biodetection, samples="SC", plottype="persample", toplot="protein_coding")
explo.plot(biodetection, samples="DUO", plottype="persample", toplot="protein_coding")
dev.off()

# COUNT DISTRIBUTION:
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/CountsBio

#counts <- dat(howell.noiseq, factor="tissue", type="countsbio")

#pdf("RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_CountDistribution.pdf", width=10, height=6)
#explo.plot(counts, samples="TI", toplot="global", plottype="boxplot")
#explo.plot(counts, samples="SC", toplot="global", plottype="boxplot")
#explo.plot(counts, samples="AC", toplot="global", plottype="boxplot")
#dev.off()

counts <- dat(re.noiseq, factor="tissue", type="countsbio")

pdf("RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_CountDistribution.pdf", width=10, height=6)
explo.plot(counts, samples="TI", toplot="global", plottype="boxplot")
explo.plot(counts, samples="SC", toplot="global", plottype="boxplot")
dev.off()

counts <- dat(af.noiseq, factor="tissue", type="countsbio")

pdf("RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_CountDistribution.pdf", width=10, height=6)
explo.plot(counts, samples="TI", toplot="global", plottype="boxplot")
explo.plot(counts, samples="SC", toplot="global", plottype="boxplot")
explo.plot(counts, samples="DUO", toplot="global", plottype="boxplot")
dev.off()

# PROTEIN CODING COUNT DISTRIBUTION:

#h.counts <- dat(howell.noiseq, factor=NULL, type="countsbio")

#pdf("RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_ProteinCoding_CountDistribution.pdf", width=20, height=5)
#explo.plot(h.counts, toplot="protein_coding", samples=NULL, plottype="boxplot")
#dev.off()

r.counts <- dat(re.noiseq, factor=NULL, type="countsbio")

pdf("RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_ProteinCoding_CountDistribution.pdf", width=30, height=5)
explo.plot(r.counts, toplot="protein_coding", samples=NULL, plottype="boxplot")
dev.off()

a.counts <- dat(af.noiseq, factor=NULL, type="countsbio")

pdf("RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_ProteinCoding_CountDistribution.pdf", width=30, height=5)
explo.plot(a.counts, toplot="protein_coding", samples=NULL, plottype="boxplot")
dev.off()

# PERCENTAGE LOW COUNT FEATURES PER SAMPLE:

#pdf("RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_ProteinCoding_PercentLowCount.pdf", width=20, height=5)
#explo.plot(h.counts, toplot="protein_coding", samples=NULL, plottype="barplot")
#dev.off()

pdf("RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_ProteinCoding_PercentLowCount.pdf", width=30, height=5)
explo.plot(r.counts, toplot="protein_coding", samples=NULL, plottype="barplot")
dev.off()

pdf("RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_ProteinCoding_PercentLowCount.pdf", width=30, height=5)
explo.plot(a.counts, toplot="protein_coding", samples=NULL, plottype="barplot")
dev.off()

# CHECK FOR LENGTH BIAS (significant p-value + R2 > 70%):
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/lengthbias

#len <- dat(howell.noiseq, factor="tissue", type="lengthbias")

#pdf("RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_ProteinCoding_LengthBias.pdf", width=7, height=5)
#explo.plot(len, samples=NULL, toplot="protein_coding")
#dev.off()

len <- dat(re.noiseq, factor="tissue", type="lengthbias")

pdf("RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_ProteinCoding_LengthBias.pdf", width=7, height=5)
explo.plot(len, samples=NULL, toplot="protein_coding")
dev.off()

len <- dat(af.noiseq, factor="tissue", type="lengthbias")

pdf("RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_ProteinCoding_LengthBias.pdf", width=7, height=5)
explo.plot(len, samples=NULL, toplot="protein_coding")
dev.off()

# CHECK FOR GC BIAS (significant p-value + R2 > 70%):
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/GCbias

#raw.gc <- dat(howell.noiseq, factor="tissue", type="GCbias")

#pdf("RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_ProteinCoding_GCBias.pdf", width=7, height=5)
#explo.plot(raw.gc, samples=NULL, toplot="protein_coding")
#dev.off()

raw.gc <- dat(re.noiseq, factor="tissue", type="GCbias")

pdf("RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_ProteinCoding_GCBias.pdf", width=7, height=5)
explo.plot(raw.gc, samples=NULL, toplot="protein_coding")
dev.off()

raw.gc <- dat(af.noiseq, factor="tissue", type="GCbias")

pdf("RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_ProteinCoding_GCBias.pdf", width=7, height=5)
explo.plot(raw.gc, samples=NULL, toplot="protein_coding")
dev.off()


#############################
# CREATE DESEQ2 QC DATASETS #
#############################

#howell.factors <- pheno.howell
#howell.factors$group <- paste0(howell.factors$tissue,".",howell.factors$diagnosis)
#rownames(howell.factors) <- howell.factors$rnaseq.id
#howel.factors$diagnosis <- as.factor(howel.factors$diagnosis)

re.factors <- pheno.re
re.factors$group <- paste0(re.factors$tissue,".",re.factors$diagnosis)
re.factors$group2 <- paste0(re.factors$group,".",re.factors$treatment)
re.factors$group3 <- paste0(re.factors$diagnosis,".",re.factors$treatment)
rownames(re.factors) <- re.factors$rnaseq.id
re.factors$diagnosis <- as.factor(re.factors$diagnosis)

af.factors <- pheno.af
af.factors$group <- paste0(af.factors$tissue,".",af.factors$diagnosis)
af.factors$group2 <- paste0(af.factors$group,".",af.factors$treatment)
af.factors$group3 <- paste0(af.factors$diagnosis,".",af.factors$treatment)
rownames(af.factors) <- af.factors$rnaseq.id
af.factors$diagnosis <- as.factor(af.factors$diagnosis)

#howell.dds <- DESeqDataSetFromMatrix(countData=raw.howell, colData=howell.factors, design=~diagnosis)
#howell.dds$diagnosis <- relevel(howell.dds$diagnosis, ref="Control")
#howell.dds <- howell.dds[rowSums(counts(howell.dds))>1,]
#howell.vst <- vst(howell.dds,blind=F)

re.dds <- DESeqDataSetFromMatrix(countData=raw.re, colData=re.factors, design=~group2)
re.dds$diagnosis <- relevel(re.dds$diagnosis, ref="Control")
re.dds <- re.dds[rowSums(counts(re.dds))>1,]
re.vst <- vst(re.dds,blind=F)

af.dds <- DESeqDataSetFromMatrix(countData=raw.af, colData=af.factors, design=~group2)
af.dds$diagnosis <- relevel(af.dds$diagnosis, ref="Control")
af.dds <- af.dds[rowSums(counts(af.dds))>1,]
af.vst <- vst(af.dds,blind=F)


###################
# DESEQ2 QC PLOTS #
###################

# EUCLDIDEAN DISTANCE HEATMAP:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-sample-to-sample-distances

#raw.dist <- dist(t(assay(howell.vst)))

#pdf("RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_RawEuclideanDistanceHeatmap.pdf", width=25, height=25)
#matrix <- as.matrix(raw.dist)
#rownames(matrix) <- howell.vst$diagnosis
#colnames(matrix) <- howell.vst$tissue
#pheatmap(matrix, clustering_distance_rows=raw.dist, clustering_distance_cols=raw.dist, col=cols, cellwidth=10, cellheight=10)
#dev.off()

raw.dist <- dist(t(assay(re.vst)))

pdf("RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_RawEuclideanDistanceHeatmap.pdf", width=20, height=20)
matrix <- as.matrix(raw.dist)
rownames(matrix) <- re.vst$group3
colnames(matrix) <- re.vst$tissue
pheatmap(matrix, clustering_distance_rows=raw.dist, clustering_distance_cols=raw.dist, col=cols, cellwidth=10, cellheight=10)
dev.off()

raw.dist <- dist(t(assay(af.vst)))

pdf("RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_RawEuclideanDistanceHeatmap.pdf", width=20, height=20)
matrix <- as.matrix(raw.dist)
rownames(matrix) <- af.vst$group3
colnames(matrix) <- af.vst$tissue
pheatmap(matrix, clustering_distance_rows=raw.dist, clustering_distance_cols=raw.dist, col=cols, cellwidth=10, cellheight=10)
dev.off()

# COOK'S CUTOFF FOR OUTLIERS:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#approach-to-count-outliers

#h.dds <- DESeq(howell.dds)
#save(howell.dds,howell.vst,h.dds,pheno.howell,file="RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_RawDataQCInput.RData")

#pdf("RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_CooksCutoff_RawNoBatchCorrection.pdf", width=20, height=7)
#boxplot(log10(assays(h.dds)[["cooks"]]), range=0, las=2)
#dev.off()

r.dds <- DESeq(re.dds)
save(re.dds,re.vst,r.dds,pheno.re,file="RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_RawDataQCInput.RData")

pdf("RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_CooksCutoff_RawNoBatchCorrection.pdf", width=30, height=7)
boxplot(log10(assays(r.dds)[["cooks"]]), range=0, las=2)
dev.off()

a.dds <- DESeq(af.dds)
save(af.dds,af.vst,a.dds,pheno.af,file="RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_RawDataQCInput.RData")

pdf("RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_CooksCutoff_RawNoBatchCorrection.pdf", width=30, height=7)
boxplot(log10(assays(a.dds)[["cooks"]]), range=0, las=2)
dev.off()

# PCA PLOTS:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#principal-component-plot-of-the-samples

#howell <- counts(h.dds, normalized=T)

#readme.howell.norm <- data.frame(Dataset=c("howell","pheno.howell","info"),Description=c("E-MTAB-5464 (Howell) normalised counts (n=76)","Howell sample info, updated 29Sep22","Ensembl GRCh37 release 107 gene annotation info"))
#save(howell,pheno.howell,info,readme.howell.norm,file="NORMALISED_COUNTS/Howell_GRCh37_NormalisedCounts.RData")

#pca.raw <- plotPCA(howell.vst, intgroup=c("group"), returnData=T)
#percentvar <- round(100*attr(pca.raw, "percentVar"))

#pdf("RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_Raw_Tissue-Diagnosis_PC1vPC2.pdf", width=10, height=10)
#ggplot(pca.raw, aes(PC1, PC2, color=group))+geom_point(size=3)+xlab(paste0("PC1: ",percentvar[1],"% variance"))+ylab(paste0("PC2: ",percentvar[2],"%variance"))
#dev.off()

#h.pca <- prcomp(t(howell))
#h.pca <- as.data.frame(h.pca$x)[,1:10]
#h.pca$rnaseq.id <- rownames(h.pca)
#h.pca <- join(pheno.howell, h.pca, type="left", match="all")
#write.table(h.pca, "NORMALISED_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_ALL_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

#pdf("NORMALISED_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_ALL_NormalisedCounts_PC1vPC2.pdf", width=15, height=5)
#p1 <- ggplot(h.pca, aes(PC1, PC2, fill=diagnosis))
#p2 <- ggplot(h.pca, aes(PC1, PC2, fill=tissue))
#p3 <- ggplot(h.pca, aes(PC1, PC2, fill=sex))
#grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis6,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite2,p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=3)
#dev.off()

re <- counts(r.dds, normalized=T)

readme.re.norm <- data.frame(Dataset=c("re","pheno.re","info"),Description=c("NGS-K.Nayak-40706 (RE) normalised counts (n=30)","RE sample info, updated 24Oct22","Ensembl GRCh37 release 107 gene annotation info"))
save(re,pheno.re,info,readme.re.norm,file="NORMALISED_COUNTS/RE_GRCh37_NormalisedCounts.RData")

pca.raw <- plotPCA(re.vst, intgroup=c("group2"), returnData=T)
percentvar <- round(100*attr(pca.raw, "percentVar"))

pdf("RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_Raw_Tissue-Diagnosis-Treatment_PC1vPC2.pdf", width=10, height=10)
ggplot(pca.raw, aes(PC1, PC2, color=group2))+geom_point(size=3)+xlab(paste0("PC1: ",percentvar[1],"% variance"))+ylab(paste0("PC2: ",percentvar[2],"%variance"))
dev.off()

r.pca <- prcomp(t(re))
r.pca <- as.data.frame(r.pca$x)[,1:10]
r.pca$rnaseq.id <- rownames(r.pca)
r.pca <- join(pheno.re, r.pca, type="left", match="all")
write.table(r.pca, "NORMALISED_COUNTS_QC/RE_GRCh37/RE_GRCh37_ALL_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/RE_GRCh37/RE_GRCh37_ALL_NormalisedCounts_PC1vPC2.pdf", width=10, height=10)
p1 <- ggplot(r.pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(r.pca, aes(PC1, PC2, fill=tissue))
p3 <- ggplot(r.pca, aes(PC1, PC2, fill=treatment))
p4 <- ggplot(r.pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite2,p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_treat,p4+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=2)
dev.off()

af <- counts(a.dds, normalized=T)

readme.af.norm <- data.frame(Dataset=c("af","pheno.af","info"),Description=c("NGS-K.Nayak-40750 (AF) normalised counts (n=90)","AF sample info, updated 24Oct22","Ensembl GRCh37 release 107 gene annotation info"))
save(af,pheno.af,info,readme.af.norm,file="NORMALISED_COUNTS/AF_GRCh37_NormalisedCounts.RData")

pca.raw <- plotPCA(af.vst, intgroup=c("group2"), returnData=T)
percentvar <- round(100*attr(pca.raw, "percentVar"))

pdf("RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_Raw_Tissue-Diagnosis-Treatment_PC1vPC2.pdf", width=10, height=10)
ggplot(pca.raw, aes(PC1, PC2, color=group2))+geom_point(size=3)+xlab(paste0("PC1: ",percentvar[1],"% variance"))+ylab(paste0("PC2: ",percentvar[2],"%variance"))
dev.off()

a.pca <- prcomp(t(af))
a.pca <- as.data.frame(a.pca$x)[,1:10]
a.pca$rnaseq.id <- rownames(a.pca)
a.pca <- join(pheno.af, a.pca, type="left", match="all")
write.table(a.pca, "NORMALISED_COUNTS_QC/AF_GRCh37/AF_GRCh37_ALL_MetaData_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("NORMALISED_COUNTS_QC/AF_GRCh37/AF_GRCh37_ALL_NormalisedCounts_PC1vPC2.pdf", width=10, height=10)
p1 <- ggplot(a.pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(a.pca, aes(PC1, PC2, fill=tissue))
p3 <- ggplot(a.pca, aes(PC1, PC2, fill=treatment))
p4 <- ggplot(a.pca, aes(PC1, PC2, fill=sex))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_sampsite2,p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_treat,p4+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("F","M")),ncol=2)
dev.off()


###########################
# CHECK FOR BATCH EFFECTS #
###########################

# Create a matrix of the normalised counts and a dataframe of the factors for sva:

#h.norm <- counts(h.dds, normalized=T)
#h.svafactors <- as.data.frame(howell.factors)

r.norm <- counts(r.dds, normalized=T)
r.svafactors <- as.data.frame(re.factors)

a.norm <- counts(a.dds, normalized=T)
a.svafactors <- as.data.frame(af.factors)

# Create models (mod=full model, mod0=null, NB, including known batch effects: diagnosis, tissue & treatment):

#h.mod <- model.matrix(~as.factor(group), data=h.svafactors)
#h.mod0 <- model.matrix(~1, data=h.svafactors)

r.mod <- model.matrix(~as.factor(group2), data=r.svafactors)
r.mod0 <- model.matrix(~1, data=r.svafactors)

a.mod <- model.matrix(~as.factor(group2), data=a.svafactors)
a.mod0 <- model.matrix(~1, data=a.svafactors)

# Calculate the first 5 SVs:

#h.svseq <- svaseq(h.norm,h.mod,h.mod0,n.sv=5)
r.svseq <- svaseq(r.norm,r.mod,r.mod0,n.sv=5)
a.svseq <- svaseq(a.norm,a.mod,a.mod0,n.sv=5)

# Plot variables against the predicted surrogate variables:

#meta_categorical <- data.frame(tissue=pheno.howell$tissue,diagnosis=pheno.howell$diagnosis,sex=pheno.howell$sex)
#meta_continuous <- data.frame(age=pheno.howell$age)
#meta_sva <- data.frame(SV1=h.svseq$sv[,1],SV2=h.svseq$sv[,2],SV3=h.svseq$sv[,3],SV4=h.svseq$sv[,4],SV5=h.svseq$sv[,5])

#h.cor <- SVcor(meta_categorical,meta_continuous,meta_sva)

#pdf("RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_SVACorrelations_5SVs.pdf", width=7, height=5)
#h.cor$plot.input
#dev.off()

#write.table(h.cor$sv.correlation, "RAW_COUNTS_QC/HOWELL_GRCh37/Howell_GRCh37_SVACorrelations_5SVs.txt", row.names=F, sep="\t", quote=F)

meta_categorical <- data.frame(tissue=pheno.re$tissue,diagnosis=pheno.re$diagnosis,sex=pheno.re$sex,treatment=pheno.re$treatment)
meta_continuous <- data.frame(age=pheno.re$age,passage=as.numeric(gsub("P","",pheno.re$passage)))
meta_sva <- data.frame(SV1=r.svseq$sv[,1],SV2=r.svseq$sv[,2],SV3=r.svseq$sv[,3],SV4=r.svseq$sv[,4],SV5=r.svseq$sv[,5])

r.cor <- SVcor(meta_categorical,meta_continuous,meta_sva)

pdf("RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_SVACorrelations_5SVs.pdf", width=8, height=7)
r.cor$plot.input
dev.off()

write.table(r.cor$sv.correlation, "RAW_COUNTS_QC/RE_GRCh37/RE_GRCh37_SVACorrelations_5SVs.txt", row.names=F, sep="\t", quote=F)

meta_categorical <- data.frame(tissue=pheno.af$tissue,diagnosis=pheno.af$diagnosis,sex=pheno.af$sex,treatment=pheno.af$treatment)
meta_continuous <- data.frame(age=pheno.af$age)
meta_sva <- data.frame(SV1=a.svseq$sv[,1],SV2=a.svseq$sv[,2],SV3=a.svseq$sv[,3],SV4=a.svseq$sv[,4],SV5=a.svseq$sv[,5])

a.cor <- SVcor(meta_categorical,meta_continuous,meta_sva)

pdf("RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_SVACorrelations_5SVs.pdf", width=8, height=7)
a.cor$plot.input
dev.off()

write.table(a.cor$sv.correlation, "RAW_COUNTS_QC/AF_GRCh37/AF_GRCh37_SVACorrelations_5SVs.txt", row.names=F, sep="\t", quote=F)


################
# SESSION INFO #
################

sessionInfo()
rm(list=ls())

