# QC OF VANCAMELBEKE GSE75214 EXPRESSION ARRAY DATA
# Background correction, normalisation, log transformation & arrayQualityMetricsQC
# Affy hugene 1.0 ST array
# 07Feb22 (fp215)
#
# Folders required:
#	- 01_QC_PLOTS
#	- OUTPUT_R_DATASETS
# 
# Files required:
#	- /home/fp215/rds/rds-fp215-working/VANCAMELBEKE/GSE75214_CEL (folder containing all raw CEL files)
#	- /home/fp215/rds/rds-fp215-working/VANCAMELBEKE/GSE75214_SampleInfo.txt (available sample info)
# 
# Useful documentation:
#	- arrayQualityMetrics: https://bioconductor.org/packages/release/bioc/vignettes/arrayQualityMetrics/inst/doc/arrayQualityMetrics.pdf
# 	- affy: https://www.bioconductor.org/packages/release/bioc/html/affy.html
#	- affycoretools: https://www.bioconductor.org/packages/release/bioc/vignettes/affycoretools/inst/doc/RefactoredAffycoretools.html
#	- annotate: https://www.bioconductor.org/packages/release/bioc/html/annotate.html
#	- pd.hugene.1.0.st: https://bioconductor.org/packages/release/data/annotation/html/pd.hugene.1.0.st.v1.html
#
# Sources useful plotting functions: PlotTools.v2.0.R 


#######################
# WORKING ENVIRONMENT #
#######################

setwd("/home/fp215/rds/rds-fp215-working/VANCAMELBEKE/GSE75214_CEL")
library(arrayQualityMetrics)
library(affy)
library(affycoretools)
library(annotate)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
library(plyr)
library(lattice)
library(ggplot2)
library(reshape)
source("/home/fp215/rds/rds-fp215-working/VANCAMELBEKE/SCRIPTS/GSE75214_PlotTools.R")


#############
# LOAD DATA #
#############

# Create a GenomeFeatureSet from the raw CEL files:

celfiles <- list.celfiles()
rawdata <- read.celfiles(celfiles)

dim(rawdata)
validObject(rawdata)

info <- read.table("/home/fp215/rds/rds-fp215-working/VANCAMELBEKE/GSE75214_SampleInfo.txt", header=T, sep="\t", stringsAsFactors=F)
pheno <- data.frame(Assay.ID=paste0(info$GEO.ID,"_",info$Sample.ID,".CEL"),Sample.ID=info$Sample.ID,GEO.ID=info$GEO.ID,Tissue=info$Tissue,Diagnosis=info$Disease,Disease.Activity=info$Disease.Activity,Inflammation=info$Inflammation)

table(pheno$Tissue,pheno$Diagnosis)
table(pheno$Diagnosis,pheno$Disease.Activity)
table(pheno$Diagnosis,pheno$Inflammation)

# Extract relative log expression data without normalising for QC:

rledata <- rma(rawdata, normalize=F)


##################
# NORMALISE DATA #
##################

# Using robust multi-array average expression measure: https://www.rdocumentation.org/packages/affy/versions/1.50.0/topics/rma
# Removes background noise, normalises & log transforms the data
# Data is output in ExpressionSet format

setwd("/home/fp215/rds/rds-fp215-working/VANCAMELBEKE/")

normdata <- rma(rawdata)

dim(normdata)
validObject(normdata)


####################
# UPDATE THE PDATA #
####################

pData(rawdata)$Assay.ID <- colnames(rawdata)
pData(rawdata)$index <- NULL
pdata <- as.data.frame(pData(rawdata))
pdata <- join(pdata, pheno, type="left", match="all")
rownames(pdata) <- pdata$Assay.ID
pdata <- AnnotatedDataFrame(pdata)
phenoData(rawdata) <- pdata

pData(normdata)$Assay.ID <- colnames(normdata)
pData(normdata)$index <- NULL
pdata <- as.data.frame(pData(normdata))
pdata <- join(pdata, pheno, type="left", match="all")
rownames(pdata) <- pdata$Assay.ID
pdata <- AnnotatedDataFrame(pdata)
phenoData(normdata) <- pdata


#################
# SAVE DATASETS #
#################

pdata <- pData(normdata)

readme.raw <- data.frame(Dataset=c("rawdata","rledata","normdata","pdata","info"),Description=c("GeneFeatureSet of raw Affy Human Gene ST 1.0 array data from colon & ileum, n=194 (Vancamelbeke, GSE75214)","ExpressionSet of RLE transformed rawdata (no normalisation)","ExpressionSet of RMA background corrected, normalised & transformed rawdata","Updated sample information table, including previously unmatched samples","Full sample information available from GSE75214"))
save(rawdata,rledata,normdata,pdata,info,readme.raw,file="Vancamelbeke_GSE75214_Raw_RMANorm_ExpressionSets.RData")

readme.raw


###########################
# QC: ARRAYQUALITYMETRICS #
###########################

# Run arrayQualityMetrics on raw & normalised datasets:

arrayQualityMetrics(expressionset=rawdata, outdir="01_QC_PLOTS/VancamelBeke_GSE75214_Raw_QC_Report", force=T, intgroup=c("Tissue","Diagnosis","Disease.Activity"), do.logtransform=T)
warnings()
arrayQualityMetrics(expressionset=normdata, outdir="01_QC_PLOTS/VancamelBeke_GSE75214_RMANorm_QC_Report", force=T, intgroup=c("Tissue","Diagnosis","Disease.Activity"), do.logtransform=T)
warnings()


###########
# QC: PCA #
###########

# PCA of RMA normalised data:
# https://www.rdocumentation.org/packages/stats/versions/3.6.1/topics/prcomp

pca <- prcomp(t(exprs(normdata)), rank=10)
pca <- as.data.frame(pca$x)
pca$Array.ID <- rownames(pca)
pca <- join(pca, pData(normdata), type="left", match="all")
write.table(pca, "01_QC_PLOTS/VancamelBeke_GSE75214_RMANorm_PC.txt", row.names=F, sep="\t", quote=F)

pdf("01_QC_PLOTS/VancamelBeke_GSE75214_RMANorm_LatticePCA_Diagnosis_n194.pdf")
splom(~pca[c("PC1","PC2","PC3","PC4","PC5")], groups=pca$Diagnosis, col=c("skyblue","grey75","yellow"), pch=rep(19,3), pscales=0, key=list(space="bottom", points=list(pch=rep(19,3), col=c("skyblue","grey75","yellow")), text=list(c("CD","Control","UC"))))
dev.off()

pdf("01_QC_PLOTS/VancamelBeke_GSE75214_RMANorm_LatticePCA_Tissue_n194.pdf")
splom(~pca[c("PC1","PC2","PC3","PC4","PC5")], groups=pca$Tissue, col=c("darkolivegreen","purple4"), pch=rep(19,2), pscales=0, key=list(space="bottom", points=list(pch=rep(19,2), col=c("darkolivegreen","purple4")), text=list(c("Colon","Ileum"))))
dev.off()

pdf("01_QC_PLOTS/VancamelBeke_GSE75214_RMANorm_LatticePCA_DiseaseActivity_n194.pdf")
splom(~pca[c("PC1","PC2","PC3","PC4","PC5")], groups=pca$Disease.Activity, col=c("black","grey50","grey75"), pch=rep(19,3), pscales=0, key=list(space="bottom", points=list(pch=rep(19,3), col=c("black","grey50","grey75")), text=list(c("active","inactive","normal"))))
dev.off()

pdf("01_QC_PLOTS/VancamelBeke_GSE75214_RMANorm_LatticePCA_Inflammation_n194.pdf")
splom(~pca[c("PC1","PC2","PC3","PC4","PC5")], groups=pca$Inflammation, col=c("grey75","darkred"), pch=rep(19,2), pscales=0, key=list(space="bottom", points=list(pch=rep(19,2), col=c("grey75","darkred")), text=list(c("No inflammation","Inflammation"))))
dev.off()


#######################
# QC: HEAT SCREE PLOT #
#######################

# Heat scree scaled PCA plot function fields: Loadings, Importance, right_marg, left_marg
# Using self-defined scaling of PCA results rather than the default

PCs_to_view <- 10

pca <- prcomp(t(exprs(normdata)))
Loadings <- as.data.frame(pca$x)
vars <- pca$sdev^2
Importance <- vars/sum(vars)

pdata <- pData(normdata)
pdata$IBD <- ifelse(pdata$Diagnosis=="control","control","IBD")
pdata$CD <- ifelse(pdata$Diagnosis=="CD","CD","control")
pdata$CD[pdata$Diagnosis=="UC"] <- NA
pdata$UC <- ifelse(pdata$Diagnosis=="UC","UC","control")
pdata$UC[pdata$Diagnosis=="CD"] <- NA
pdata$Inflammation.2 <- ifelse(pdata$Inflammation=="Y",1,0)
pdata$Activity.2 <- ifelse(pdata$Disease.Activity=="active",1,0)
pdata$Activity.2[pdata$Diagnosis=="control"] <- NA

meta_categorical <- pdata[c("Tissue","IBD","CD","UC","Inflammation.2","Activity.2")]
colnames(meta_categorical) <- c("Tissue","IBD/Control","CD/Control","UC/Control","Inflammation (all)","Activity (IBD)")
ord <- 1:length(c(colnames(meta_categorical)))

pdf("01_QC_PLOTS/VancamelBeke_GSE75214_RMANorm_HeatScreePlot_n194.pdf", width=9, height=6)
heat.scree.cat(Loadings, Importance, 2.5, 2.7)
dev.off()


###############################
# QC: RELATIVE LOG EXPRESSION #
###############################

library(tidy)

# Calculate medians:

meds <- rowMedians(as.matrix(exprs(rledata)))
rle.mat <- sweep(exprs(rledata), 1, meds)
rle.mat <- as.data.frame(rle.mat)
rle.gath <- gather(rle.mat, patient_array, log2_expression_deviation)

# Check the limits of the data:

summary(rle.gath$log2_expression_deviation)

# Create boxplot:

pdf("01_QC_PLOTS/VancamelBeke_GSE75214_Raw_RLEBoxplot_n194.pdf", width=30, height=5)
ggplot(rle.gath, aes(patient_array, log2_expression_deviation))+geom_boxplot(outlier.shape=NA)+ylim(c(-9,9))+theme(axis.text.x=element_text(angle=60,size=6.5,hjust=1,face="bold"))+geom_hline(yintercept=0, colour="blue")+geom_hline(yintercept=-2, colour="red", linetype=3)+geom_hline(yintercept=2, colour="red", linetype=3)
dev.off()


################
# SESSION INFO #
################

sessionInfo()
rm(list=ls())

