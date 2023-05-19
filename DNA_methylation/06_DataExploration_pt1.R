# DATA EXPLORATION (PART 1):
# Folders needed: PCA_PLOTS, HEAT_SCREE_PLOTS, CLUSTER_PLOTS
# 10Oct22


# SET UP ENVIRONMENT & LOAD DATA:

library(plyr)
library(ggplot2)
library(reshape)
library(scales)
library(ClassDiscovery)
library(rafalib)
library(dendextend)
require(gridExtra)
library(grid)
source("/home/fp215/SCRIPTS/ADDENBROOKES_SCRIPTS/GENERAL_PURPOSE/PlotTools.v3.4.R")

PCs_to_view <- 10

# Dataset:

load("AllEpicOrganoid_WorkingDataset_10Oct22.RData")

readme.epic

# Load EPIC array batch & passage colour schemes:

load("ArrayColourSchemes_65batches.RData")
load("PassageColourSchemes_P1-P8_P10_P16.RData")

# Plot themes to match Rachel Edgar's methylation score plots:

myTheme_meanmeth_label <- theme(legend.position="bottom", axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_blank(), panel.border=element_blank(), legend.text=element_text(size=rel(2)))
myTheme_meanmeth_meth <- theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())


# CLINICAL DATA ANNOTATION UPDATES:

cont <- c("T474","T477","T493","T528","T640")
samples$DUO.Inflammation[samples$case.no %in% cont] <- "N"
samples$TI.Inflammation[samples$case.no %in% cont] <- "N"
samples$SC.Inflammation[samples$case.no %in% cont] <- "N"
samples$Biologics[samples$case.no %in% cont] <- "N"
samples$Surgery[samples$case.no %in% cont] <- "N"
samples$AZA[samples$case.no %in% cont] <- "N"
samples$Perianal.Disease[samples$case.no %in% cont] <- "N"
samples$Disease.Severity[samples$case.no %in% cont] <- "control"

save(samples,beta,mval,mhc,beta.mhc,beta.nlrc5,gene.info,readme.epic,file="AllEpicOrganoid_WorkingDataset_10Oct22.RData")


# ADD NUMERIC ANNOTATION FOR CLINICAL MEASUREMENTS:
# TI.Inflammation/SC.Inflammation/DUO.Inflammation: 1/0 = yes (macroscopic/microscopic)/ no
# Biologics/Surgery/AZA/Perianal Disease: 1/0 = yes/ no
# Disease Severity: 0/1/2/3 = control/ mild/ moderate/ severe
# NB, T490 is currently missing TI inflammation status - list as "N" for now, until Federica can update
# NB2, T338 is currently missing perianal disease status - list as "N" for now, until Federica can update

samples$TI.numeric <- ifelse(samples$TI.Inflammation=="N",0,1)
samples$TI.numeric[samples$case.no=="T490"] <- 0
samples$SC.numeric <- ifelse(samples$SC.Inflammation=="N",0,1)
samples$DUO.numeric <- ifelse(samples$DUO.Inflammation=="N",0,1)
samples$Biologics.numeric <- ifelse(samples$Biologics=="Y",1,0)
samples$Biologics.numeric[samples$diagnosis=="Control"] <- 0
samples$Surgery.numeric <- ifelse(samples$Surgery=="Y",1,0)
samples$Surgery.numeric[samples$diagnosis=="Control"] <- 0
samples$AZA.numeric <- ifelse(samples$AZA=="Y",1,0)
samples$AZA.numeric[samples$diagnosis=="Control"] <- 0
samples$Perianal.numeric <- ifelse(samples$Perianal.Disease=="Y",1,0)
samples$Perianal.numeric[samples$case.no=="T338"] <- 0
samples$Perianal.numeric[samples$diagnosis=="Control"] <- 0
samples$Severity <- samples$Disease.Severity
samples$Severity <- gsub("control",0,samples$Severity)
samples$Severity <- gsub("mild",1,samples$Severity)
samples$Severity <- gsub("moderate",2,samples$Severity)
samples$Severity <- gsub("severe",3,samples$Severity)


# ADD COLOUR SCHEME:

samples$col.diagnosis <- as.factor(samples$diagnosis)
levels(samples$col.diagnosis) <- c("dodgerblue3","lightgrey","#636B4D","darkgoldenrod1")
samples$col.diagnosis <- as.character(samples$col.diagnosis)


# SPLIT DATASET BY TISSUE:

ti.samples <- subset(samples, tissue=="TI")
ti.beta <- beta[,ti.samples$array.id]
ti.mval <- mval[,ti.samples$array.id]
ti.mhc <- beta.mhc[,ti.samples$array.id]

sc.samples <- subset(samples, tissue=="SC")
sc.beta <- beta[,sc.samples$array.id]
sc.mval <- mval[,sc.samples$array.id]
sc.mhc <- beta.mhc[,sc.samples$array.id]

duo.samples <- subset(samples, tissue=="DUO")
duo.beta <- beta[,duo.samples$array.id]
duo.mval <- mval[,duo.samples$array.id]
duo.mhc <- beta.mhc[,duo.samples$array.id]

nrow(ti.samples)
dim(ti.beta)
dim(ti.mval)
dim(ti.mhc)
nrow(sc.samples)
dim(sc.beta)
dim(sc.mval)
dim(sc.mhc)
nrow(duo.samples)
dim(duo.beta)
dim(duo.mval)
dim(duo.mhc)

save(ti.samples,ti.beta,ti.mval,ti.mhc,gene.info,file="AllEPICOrganoid_TI_WorkingDataset_10Oct22.RData")
save(sc.samples,sc.beta,sc.mval,sc.mhc,gene.info,file="AllEPICOrganoid_SC_WorkingDataset_10Oct22.RData")
save(duo.samples,duo.beta,duo.mval,duo.mhc,gene.info,file="AllEPICOrganoid_DUO_WorkingDataset_10Oct22.RData")


# PCA PLOTS:

passage$Passage <- as.factor(passage$numericP)

ti.pca <- prcomp(t(ti.beta))
ti.pca <- as.data.frame(ti.pca$x)[,1:10]
ti.pca$array.id <- rownames(ti.pca)
ti.pca <- join(ti.pca, ti.samples, type="left", match="all")
write.table(ti.pca, "PCA_PLOTS/PaediatricOrganoid_TI_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("PCA_PLOTS/PaediatricOrganoid_TI_PC1vPC2.pdf", width=10, height=10)
p1 <- ggplot(ti.pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(ti.pca, aes(PC1, PC2, fill=as.factor(passage)))
p3 <- ggplot(ti.pca, aes(PC1, PC2, fill=sex))
p4 <- ggplot(ti.pca, aes(PC1, PC2, fill=RE.batch))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_RE_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Passage",values=passage$cols,labels=passage$Passage),p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("Female","Male")),p4+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Batch",values=c("green4","goldenrod","darkorange4"),labels=c("NEW","OLD","RE")),ncol=2)
dev.off()

sc.pca <- prcomp(t(sc.beta))
sc.pca <- as.data.frame(sc.pca$x)[,1:10]
sc.pca$array.id <- rownames(sc.pca)
sc.pca <- join(sc.pca, sc.samples, type="left", match="all")
write.table(sc.pca, "PCA_PLOTS/PaediatricOrganoid_SC_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("PCA_PLOTS/PaediatricOrganoid_SC_PC1vPC2.pdf", width=10, height=10)
p1 <- ggplot(sc.pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(sc.pca, aes(PC1, PC2, fill=as.factor(passage)))
p3 <- ggplot(sc.pca, aes(PC1, PC2, fill=sex))
p4 <- ggplot(sc.pca, aes(PC1, PC2, fill=RE.batch))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_RE_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Passage",values=passage$cols,labels=passage$Passage),p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("Female","Male")),p4+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Batch",values=c("green4","goldenrod","darkorange4"),labels=c("NEW","OLD","RE")),ncol=2)
dev.off()

duo.pca <- prcomp(t(duo.beta))
duo.pca <- as.data.frame(duo.pca$x)[,1:10]
duo.pca$array.id <- rownames(duo.pca)
duo.pca <- join(duo.pca, duo.samples, type="left", match="all")
write.table(duo.pca, "PCA_PLOTS/PaediatricOrganoid_DUO_PCs.txt", row.names=F, sep="\t", quote=F)

pdf("PCA_PLOTS/PaediatricOrganoid_DUO_PC1vPC2.pdf", width=10, height=10)
p1 <- ggplot(duo.pca, aes(PC1, PC2, fill=diagnosis))
p2 <- ggplot(duo.pca, aes(PC1, PC2, fill=as.factor(passage)))
p3 <- ggplot(duo.pca, aes(PC1, PC2, fill=sex))
p4 <- ggplot(duo.pca, aes(PC1, PC2, fill=RE.batch))
grid.arrange(p1+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+fillscale_RE_diagnosis2,p2+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Passage",values=passage$cols,labels=passage$Passage),p3+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Sex",values=c("red","blue"),labels=c("Female","Male")),p4+geom_point(shape=21,size=2,color="black")+theme_bw()+theme(legend.title=element_text(face="bold"))+scale_fill_manual(name="Batch",values=c("green4","goldenrod","darkorange4"),labels=c("NEW","OLD","RE")),ncol=2)
dev.off()


# HEAT SCREE PLOTS:
# NB, currently, T490 has no annotation for the "TI.Inflammation" column - temporarily label as "N" & update if that alters

ti.pca <- prcomp(t(ti.beta))
Loadings <- as.data.frame(ti.pca$x)
vars <- ti.pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- ti.samples[c("diagnosis","sex","sentrix.id","RE.batch","Wnt.type")]
meta_continuous <- ti.samples[c("age","numeric.passage")]
colnames(meta_categorical) <- c("Diagnosis","Sex","Sentrix ID","Batch","Wnt Type")
colnames(meta_continuous) <- c("Age","Passage")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_TI_SampleInfo_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

meta_categorical <- ti.samples[c("sex","TI.numeric","Biologics.numeric","Surgery.numeric","AZA.numeric","Perianal.numeric")]
meta_continuous <- ti.samples[c("age","Severity")]
colnames(meta_categorical) <- c("Sex","Inflammation","Biologics","Surgery","AZA","Perianal Disease")
colnames(meta_continuous) <- c("Age","Disease Severity")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_TI_OutcomeMeasures_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

sc.pca <- prcomp(t(sc.beta))
Loadings <- as.data.frame(sc.pca$x)
vars <- sc.pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- sc.samples[c("diagnosis","sex","sentrix.id","RE.batch","Wnt.type")]
meta_continuous <- sc.samples[c("age","numeric.passage")]
colnames(meta_categorical) <- c("Diagnosis","Sex","Sentrix ID","Batch","Wnt Type")
colnames(meta_continuous) <- c("Age","Passage")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_SC_SampleInfo_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

meta_categorical <- sc.samples[c("sex","SC.numeric","Biologics.numeric","Surgery.numeric","AZA.numeric","Perianal.numeric")]
meta_continuous <- sc.samples[c("age","Severity")]
colnames(meta_categorical) <- c("Sex","Inflammation","Biologics","Surgery","AZA","Perianal Disease")
colnames(meta_continuous) <- c("Age","Disease Severity")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_SC_OutcomeMeasures_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

duo.pca <- prcomp(t(duo.beta))
Loadings <- as.data.frame(duo.pca$x)
vars <- duo.pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- duo.samples[c("diagnosis","sex","sentrix.id","RE.batch","Wnt.type")]
meta_continuous <- duo.samples[c("age","numeric.passage")]
colnames(meta_categorical) <- c("Diagnosis","Sex","Sentrix ID","Batch","Wnt Type")
colnames(meta_continuous) <- c("Age","Passage")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_DUO_SampleInfo_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

meta_categorical <- duo.samples[c("sex","DUO.numeric","Biologics.numeric","Surgery.numeric","AZA.numeric","Perianal.numeric")]
meta_continuous <- duo.samples[c("age","Severity")]
colnames(meta_categorical) <- c("Sex","Inflammation","Biologics","Surgery","AZA","Perianal Disease")
colnames(meta_continuous) <- c("Age","Disease Severity")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_DUO_OutcomeMeasures_HeatScree.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()


# BASIC CLUSTER PLOTS:
# NB, save input for easy redrawing

d.ti <- dist(t(ti.mhc))
hc.ti <- hclust(d.ti, method="complete")
dend.ti <- as.dendrogram(hc.ti)
order.ti <- ti.samples[match(labels(dend.ti),ti.samples$array.id),]
order.ti$order <- c(1:nrow(order.ti))
save(d.ti,hc.ti,dend.ti,order.ti,file="CLUSTER_PLOTS/PaediatricOrganoid_TI_ClusteringInput.RData")

meth <- ggplot(order.ti, aes(reorder(array.id,order), Avg.MHCI.CpG, group=1))+geom_point()+geom_line()+theme_bw()+myTheme_meanmeth_meth
vp1 <- viewport(height=unit(0.15,"npc"), width=unit(0.87,"npc"), just=c("center","top"), y=0.17, x=0.495)

pdf("CLUSTER_PLOTS/PaediatricOrganoid_TI_BasicAverageMethylationClustering.pdf", width=30, height=15, onefile=F)
par(mfcol=c(1,1), mar=c(15,5,0,5)+0.1, oma=c(0,1,0,1))
plot.new()
myplclust(hc.ti, labels=ti.samples$case.no, lab.col=ti.samples$col.diagnosis, cex=1.2, main=NA, ylab=NA)
print(meth, vp=vp1)
grid.text("Average\nMHCI\nMethylation", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.1,"npc"), gp=gpar(fontsize=15))
dev.off()

d.sc <- dist(t(sc.mhc))
hc.sc <- hclust(d.sc, method="complete")
dend.sc <- as.dendrogram(hc.sc)
order.sc <- sc.samples[match(labels(dend.sc),sc.samples$array.id),]
order.sc$order <- c(1:nrow(order.sc))
save(d.sc,hc.sc,dend.sc,order.sc,file="CLUSTER_PLOTS/PaediatricOrganoid_SC_ClusteringInput.RData")

meth <- ggplot(order.sc, aes(reorder(array.id,order), Avg.MHCI.CpG, group=1))+geom_point()+geom_line()+theme_bw()+myTheme_meanmeth_meth
vp1 <- viewport(height=unit(0.15,"npc"), width=unit(0.87,"npc"), just=c("center","top"), y=0.17, x=0.495)

pdf("CLUSTER_PLOTS/PaediatricOrganoid_SC_BasicAverageMethylationClustering.pdf", width=30, height=15, onefile=F)
par(mfcol=c(1,1), mar=c(15,5,0,5)+0.1, oma=c(0,1,0,1))
plot.new()
myplclust(hc.sc, labels=sc.samples$case.no, lab.col=sc.samples$col.diagnosis, cex=1.2, main=NA, ylab=NA)
print(meth, vp=vp1)
grid.text("Average\nMHCI\nMethylation", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.1,"npc"), gp=gpar(fontsize=15))
dev.off()

d.duo <- dist(t(duo.mhc))
hc.duo <- hclust(d.duo, method="complete")
dend.duo <- as.dendrogram(hc.duo)
order.duo <- duo.samples[match(labels(dend.duo),duo.samples$array.id),]
order.duo$order <- c(1:nrow(order.duo))
save(d.duo,hc.duo,dend.duo,order.duo,file="CLUSTER_PLOTS/PaediatricOrganoid_DUO_ClusteringInput.RData")

meth <- ggplot(order.duo, aes(reorder(array.id,order), Avg.MHCI.CpG, group=1))+geom_point()+geom_line()+theme_bw()+myTheme_meanmeth_meth
vp1 <- viewport(height=unit(0.15,"npc"), width=unit(0.85,"npc"), just=c("center","top"), y=0.2, x=0.495)

pdf("CLUSTER_PLOTS/PaediatricOrganoid_DUO_BasicAverageMethylationClustering.pdf", width=20, height=10, onefile=F)
par(mfcol=c(1,1), mar=c(15,5,0,5)+0.1, oma=c(0,1,0,1))
plot.new()
myplclust(hc.duo, labels=duo.samples$case.no, lab.col=duo.samples$col.diagnosis, cex=1.2, main=NA, ylab=NA)
print(meth, vp=vp1)
grid.text("Average\nMHCI\nMethylation", just=c("right","center"), x=unit(0.06,"npc"), y=unit(0.13,"npc"), gp=gpar(fontsize=15))
dev.off()


# SESSION INFO:

sessionInfo()


