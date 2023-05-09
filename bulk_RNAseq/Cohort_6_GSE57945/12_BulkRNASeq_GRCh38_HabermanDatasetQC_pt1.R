# QC OF HABERMAN GSE57945 RAW BULK RNASEQ COUNTS - Part 1
# NOISeq & DESeq2 QC plots, check for batch effects
# 17Jan22 (fp215)
# 
# Folders required:
#     - RAW_COUNTS_QC
#       
# Files required:
#     - Raw counts: /Users/fp215/Documents/DATASETS/RNASEQ/Haberman_GSE57945/GSE57945_COUNTS/Haberman_GSE57945_per-gene-counts.2.txt
#
# Useful documentation:
#     - biomaRt: https://bioconductor.org/packages/release/bioc/vignettes/biomaRt/inst/doc/biomaRt.html
#     - NOISeq documentation: https://bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf
#     - DESeq2 documentation: https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
#     - SVA documentation: https://bioconductor.org/packages/release/bioc/vignettes/sva/inst/doc/sva.pdf


#######################
# WORKING ENVIRONMENT #
#######################

setwd("/Users/fp215/Documents/DATASETS/RNASEQ/Haberman_GSE57945")
library(biomaRt)
library(NOISeq)
library(DESeq2)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
library(sva)

cols <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)


##############################
# EXTRACT ENSEMBL ANNOTATION #
##############################

# Extract annotation for GRCh38 version 99:

ensembl <- useEnsembl(biomart="ensembl",version=99,dataset="hsapiens_gene_ensembl")

# Pull out information of interest:

chrs <- c(1:22,"X","Y")
info <- getBM(attributes=c("external_gene_name", "ensembl_gene_id", "chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype", "transcript_count"), mart=ensembl)
info$length <- info$end_position - info$start_position
info <- rename(info, c("external_gene_name"="gene.symbol","ensembl_gene_id"="gene","chromosome_name"="chr","start_position"="start","end_position"="end","percentage_gene_gc_content"="percent.gc","gene_biotype"="biotype","transcript_count"="tc"))
info <- subset(info, chr %in% chrs)

# Save annotation info:

save(info,file="~/Documents/DATASETS/RNASEQ/REFERENCE/EnsemblGRCh38v99_Info.RData")


#############
# LOAD DATA #
#############

raw <- read.table("~/Documents/DATASETS/RNASEQ/Haberman_GSE57945/GSE57945_COUNTS/Haberman_GSE57945_per-gene-counts.2.txt", header=T, sep="\t", stringsAsFactors=F)
rownames(raw) <- raw$Geneid
raw$Geneid <- NULL

# Pull phenotype table from RPKM/TPM dataset previously produced (contains RPKM & TPM values + a phenotype table ("pheno") & a table of the samples missing from the TPM dataset:

load("/Users/fp215/Documents/PAST_PROJECT_USEFUL_STUFF/REPLECATION_DATASETS/HABERMAN_BULK_RNASEQ/GSE57945.RData")

# Match order of phenotype table to order of raw counts table:

pheno <- pheno[match(colnames(raw), pheno$SRA.ID),]

table(pheno$Gender)
table(pheno$Diagnosis)
table(pheno$Diagnosis.2)
table(pheno$Deep.Ulcer)
table(pheno$Inflammation)

# Save raw dataset:

readme.raw <- data.frame(Dataset=c("raw","pheno"),Description=c("Haberman GSE57945 raw counts (n=322)","Haberman GSE57945 sample info (n=322)"))
save(raw,pheno,readme.raw,file="GSE57945_RawCounts.RData")


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

raw.noiseq <- readData(data=raw, length=length, gc=gc, biotype=biotype, chromosome=coords, factors=pheno)


###################
# NOISEQ QC PLOTS #
###################

# BIOTYPE DISTRIBUTION:
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/Biodetection

biodetection <- dat(raw.noiseq, k=0, type="biodetection", factor="Diagnosis")

pdf("RAW_COUNTS_QC/GSE57945_BiotypeDistribution.pdf", width=10, height=5)
explo.plot(biodetection, samples="CD", plottype="persample", toplot="protein_coding")
explo.plot(biodetection, samples="UC", plottype="persample", toplot="protein_coding")
explo.plot(biodetection, samples="Control", plottype="persample", toplot="protein_coding")
dev.off()

# COUNT DISTRIBUTION:
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/CountsBio

counts <- dat(raw.noiseq, factor="Diagnosis", type="countsbio")

pdf("RAW_COUNTS_QC/GSE57945_CountDistribution.pdf", width=10, height=5)
explo.plot(counts, samples="CD", toplot="global", plottype="boxplot")
explo.plot(counts, samples="UC", toplot="global", plottype="boxplot")
explo.plot(counts, samples="Control", toplot="global", plottype="boxplot")
dev.off()

# SATURATION:
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/Saturation

saturation <- dat(raw.noiseq, k=0, ndepth=6, type="saturation")

pdf("RAW_COUNTS_QC/GSE57945_DetectedFeatures_vs_SeqDepth.pdf", width=10, height=10)
explo.plot(saturation, toplot="global", samples=1:length(pheno))
explo.plot(saturation, toplot="protein_coding", samples=1:length(pheno))
dev.off()

# PROTEIN CODING COUNT DISTRIBUTION:

counts <- dat(raw.noiseq, factor=NULL, type="countsbio")

pdf("RAW_COUNTS_QC/GSE57945_ProteinCoding_CountDistribution.pdf", width=50, height=5)
explo.plot(counts, toplot="protein_coding", samples=NULL, plottype="boxplot")
dev.off()

# PERCENTAGE LOW COUNT FEATURES PER SAMPLE:

pdf("RAW_COUNTS_QC/GSE57945_ProteinCoding_PercentLowCount.pdf", width=50, height=5)
explo.plot(counts, toplot="protein_coding", samples=NULL, plottype="barplot")
dev.off()

# CHECK FOR LENGTH BIAS (significant p-value + R2 > 70%):
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/lengthbias

len <- dat(raw.noiseq, factor="Diagnosis", type="lengthbias")

pdf("RAW_COUNTS_QC/GSE57945_ProteinCoding_LenghBias.pdf", width=7, height=5)
explo.plot(len, samples=NULL, toplot="protein_coding")
dev.off()

# CHECK FOR GC BIAS (significant p-value + R2 > 70%):
# https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/GCbias

raw.gc <- dat(raw.noiseq, factor="Diagnosis", type="GCbias")

pdf("RAW_COUNTS_QC/GSE57945_ProteinCoding_GCBias.pdf", width=7, height=5)
explo.plot(raw.gc, samples=NULL, toplot="protein_coding")
dev.off()


#############################
# CREATE DESEQ2 QC DATASETS #
#############################

factors <- pheno
rownames(factors) <- factors$SRA.ID
factors$SRA.ID <- NULL
factors <- as.matrix(factors)

raw.deseq <- DESeqDataSetFromMatrix(countData=raw, colData=factors, design=~Diagnosis)
raw.deseq$Diagnosis <- relevel(raw.deseq$Diagnosis, ref="Control")
raw.deseq <- raw.deseq[rowSums(counts(raw.deseq))>1,]

raw.deseq.vst <- vst(raw.deseq,blind=F)


###################
# DESEQ2 QC PLOTS #
###################

# EUCLDIDEAN DISTANCE HEATMAP:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#heatmap-of-the-sample-to-sample-distances

cols <- colorRampPalette(rev(brewer.pal(9,"Blues")))(255)
raw.dist <- dist(t(assay(raw.deseq.vst)))

pdf("RAW_COUNTS_QC/GSE57945_EuclideanDistanceHeatmap.pdf", width=25, height=25)
matrix <- as.matrix(raw.dist)
rownames(matrix) <- raw.deseq.vst$Diagnosis
colnames(matrix) <- NULL
pheatmap(matrix, clustering_distance_rows=raw.dist, clustering_distance_cols=raw.dist, col=cols, cellwidth=5, cellheight=5)
dev.off()

# PCA PLOTS:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#principal-component-plot-of-the-samples

pca.raw <- plotPCA(raw.deseq.vst, intgroup=c("Diagnosis"), returnData=T)
percentvar <- round(100*attr(pca.raw, "percentVar"))

pdf("RAW_COUNTS_QC/GSE57945_Diagnosis_PC1vPC2.pdf", width=10, height=10)
ggplot(pca.raw, aes(PC1, PC2, color=Diagnosis))+geom_point(size=3)+xlab(paste0("PC1: ",percentvar[1],"% variance"))+ylab(paste0("PC2: ",percentvar[2],"%variance"))
dev.off()

# COOK'S CUTOFF FOR OUTLIERS:
# https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#approach-to-count-outliers

dds <- DESeq(raw.deseq)
save(dds,file="RAW_COUNTS_QC/GSE57945_Raw_dds.RData")

pdf("RAW_COUNTS_QC/GSE57945_CooksCutoff_NoBatchCorrection.pdf", width=40, height=7)
boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2)
dev.off()


###########################
# CHECK FOR BATCH EFFECTS #
###########################

# Create a matrix of the normalised counts and a dataframe of the factors for sva:

norm <- counts(dds, normalized=T)
svafactors <- as.data.frame(factors)

# Create models (mod=full model, mod0=null, NB, including known batch effects: diagnosis & tissue):

mod <- model.matrix(~as.factor(Diagnosis), data=svafactors)
mod0 <- model.matrix(~1, data=svafactors)

# Calculate the first 2 SVs:

svseq <- svaseq(norm,mod,mod0,n.sv=2)

# Run DESeq2 correcting for 2 SVs:

dds.2 <- dds
dds.2$SV1 <- svseq$sv[,1]
dds.2$SV2 <- svseq$sv[,2]

# Create Euclidean distance plots for uncorrected & corrected data:

qc.0 <- vst(dds,blind=F)
qc.2 <- vst(dds.2,blind=F)

dist.0 <- dist(t(assay(qc.0)))
dist.2 <- dist(t(assay(qc.2)))

pdf("RAW_COUNTS_QC/GSE57945_SVA_0SVs_EuclideanDistance.pdf", width=25, height=25)
matrix <- as.matrix(dist.0)
rownames(matrix) <- qc.0$Diagnosis
colnames(matrix) <- NULL
pheatmap(matrix, clustering_distance_rows=dist.0, clustering_distance_cols=dist.0, col=cols, cellwidth=5, cellheight=5)
dev.off()

pdf("RAW_COUNTS_QC/GSE57945_SVA_2SVs_EuclideanDistance.pdf", width=25, height=25)
matrix <- as.matrix(dist.2)
rownames(matrix) <- qc.2$Diagnosis
colnames(matrix) <- NULL
pheatmap(matrix, clustering_distance_rows=dist.2, clustering_distance_cols=dist.2, col=cols, cellwidth=5, cellheight=5)
dev.off()

# Plot the factors estimated by SVA for 2 SVs:

pdf("RAW_COUNTS_QC/GSE57945_SVA_2SVs_Diagnosis_Factors.pdf")
par(mfrow = c(2, 1), mar = c(3,5,3,1))
for (i in 1:2) {
  stripchart(svseq$sv[,i] ~ dds$Diagnosis, vertical = TRUE, main = paste0("SV", i))
  abline(h = 0)
 }
dev.off()


################
# SESSION INFO #
################

sessionInfo()
