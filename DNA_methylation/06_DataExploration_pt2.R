# DATA EXPLORATION (PART 2):
# Folders needed: HEAT_SCREE_PLOTS, CLUSTER_PLOTS, MHCI_SCORE_DISTRIBUTION, MHCI_SCORE_CORRELATIONS
# Using NLRC5 paper specific plot themes & schemes in NLRC5_PlotTools.R in addition to the usual
# 31Oct22 (fp215)


##################################
# SET UP ENVIRONMENT & LOAD DATA #
##################################

library(plyr)
library(ggplot2)
library(reshape)
library(scales)
library(ClassDiscovery)
library(rafalib)
library(dendextend)
require(gridExtra)
library(grid)
library(ggpubr)
library(rstatix)
library(ggnewscale)
library(dplyr)
source("/home/fp215/SCRIPTS/ADDENBROOKES_SCRIPTS/GENERAL_PURPOSE/PlotTools.v3.4.R")
source("/home/fp215/rds/rds-fp215-working/METHYLATION_ARRAY_ANALYSIS/NLRC5_PAPER_2/SCRIPTS/NLRC5_PlotTools.R")

PCs_to_view <- 10


# DATASETS:
# Full unsplit analysis ready dataset:

load("AllEpicOrganoid_WorkingDataset_10Oct22.RData")

readme.epic

# Overwrite "samples" with updated sample info from 31Oct22 (also includes cluster groupings):
# Also update nomenclature for SC & DUO inflammation columns (TI already updated)
# Sample info updated: IBD-U samples categorised as either CD or UC (depending on if they're "CD-like" or "UC-like"):
# cdlike <- c("279","374","T013","T027","T051","T148","T153","T166","T172","T190","T338","T343","T439","T440","T454","T472","T212","T232")
# uclike <- c("266","T143")
# T454: Disease.Severity altered from "mild" to "moderate"
# T454 & T212: TI.Inflammation altered from "Normal" to "Microscopic"
# T212: TI.Inflammation altered from "Microscopic" to "Macroscopic"
# T166: DUO.Inflammation altered from "N" to "Microscopic"
# All updated to the sample info sheet externally & saved in: AllEPICOrganoid_MHCIannotated_UpdatedSampleInfo_31Oct22.txt

samples <- read.table("AllEPICOrganoid_MHCIannotated_UpdatedSampleInfo_31Oct22.txt", header=T, sep="\t", stringsAsFactors=F)
samples$SC.Inflammation <- gsub("N","Normal",samples$SC.Inflammation)
samples$DUO.Inflammation <- gsub("N","Normal",samples$DUO.Inflammation)

# Additional annotation for plotting:

samples$TI.numeric <- ifelse(samples$TI.Inflammation=="Normal",0,1)
samples$SC.numeric <- ifelse(samples$SC.Inflammation=="Normal",0,1)
samples$DUO.numeric <- ifelse(samples$DUO.Inflammation=="Normal",0,1)
samples$Biologics.numeric <- ifelse(samples$Biologics=="Y",1,0)
samples$Biologics.numeric[samples$diagnosis=="Control"] <- 0
samples$Surgery.numeric <- ifelse(samples$Surgery=="Y",1,0)
samples$Surgery.numeric[samples$diagnosis=="Control"] <- 0
samples$AZA.numeric <- ifelse(samples$AZA=="Y",1,0)
samples$AZA.numeric[samples$diagnosis=="Control"] <- 0
samples$Perianal.numeric <- ifelse(samples$Perianal.Disease=="Y",1,0)
samples$Perianal.numeric[samples$diagnosis=="Control"] <- 0
samples$Severity <- samples$Disease.Severity
samples$Severity <- gsub("control",0,samples$Severity)
samples$Severity <- gsub("mild",1,samples$Severity)
samples$Severity <- gsub("moderate",2,samples$Severity)
samples$Severity <- gsub("severe",3,samples$Severity)
samples$col.diagnosis <- as.factor(samples$diagnosis)
levels(samples$col.diagnosis) <- c("dodgerblue3","lightgrey","darkgoldenrod1")
samples$col.diagnosis <- as.character(samples$col.diagnosis)
samples$TI.MHCI.Group.2 <- samples$TI.MHCI.Group
samples$TI.MHCI.Group.2[samples$TI.Cluster==1] <- "Intermediate"
samples$TI.MHCI.Group.2[samples$TI.Cluster==2] <- "High"
samples$SC.MHCI.Group.2 <- samples$SC.MHCI.Group
samples$SC.MHCI.Group.2[samples$SC.Cluster==1] <- "High"
samples$SC.MHCI.Group.2[samples$SC.Cluster==2] <- "Intermediate"
samples$severe.v.other <- ifelse(samples$Disease.Severity=="severe","severe","other")
samples$severe.v.other[samples$diagnosis=="Control"] <- "control"
samples$mild.v.other <- ifelse(samples$Disease.Severity=="mild","mild","other")
samples$mild.v.other[samples$diagnosis=="Control"] <- "control"

# Full tissue specific analysis ready datasets:

load("AllEPICOrganoid_TI_WorkingDataset_24Oct22.RData")
load("AllEPICOrganoid_SC_WorkingDataset_24Oct22.RData")
load("AllEPICOrganoid_DUO_WorkingDataset_24Oct22.RData")

# Update the tissue specific sample info tables:

ti.samples <- subset(samples, tissue=="TI")
sc.samples <- subset(samples, tissue=="SC")
duo.samples <- subset(samples, tissue=="DUO")

all(ti.samples$array.id==colnames(ti.beta))
all(sc.samples$array.id==colnames(sc.beta))
all(duo.samples$array.id==colnames(duo.beta))

# Previously created MHCI probe clustering output:

load("CLUSTER_PLOTS/PaediatricOrganoid_TI_ClusteringInput.RData")
load("CLUSTER_PLOTS/PaediatricOrganoid_SC_ClusteringInput.RData")
load("CLUSTER_PLOTS/PaediatricOrganoid_DUO_ClusteringInput.RData")

# Load up "clust" (summarises all clustering info), update the diagnosis column & re-save:

load("CLUSTER_PLOTS/PaediatricOrganoid_ClusteringGroups.RData")
clust$diagnosis <- NULL
clust <- join(clust, unique(samples[,c("case.no","diagnosis")]), type="left", match="all")
save(clust,file="CLUSTER_PLOTS/PaediatricOrganoid_ClusteringGroups.RData")

# Re-save the datasets with updated sample info:

save(ti.samples,ti.beta,ti.mval,ti.mhc,gene.info,file="AllEPICOrganoid_TI_WorkingDataset_31Oct22.RData")
save(sc.samples,sc.beta,sc.mval,sc.mhc,gene.info,file="AllEPICOrganoid_SC_WorkingDataset_31Oct22.RData")
save(duo.samples,duo.beta,duo.mval,duo.mhc,gene.info,file="AllEPICOrganoid_DUO_WorkingDataset_31Oct22.RData")
save(samples,file="AllEPICOrganoid_MHCIannotated_UpdatedSampleInfoWithPlotInfo_31Oct22.RData")

# Additional colour schemes & updates to old colour schemes:

myColors_RE <- c("lightgrey","darkgoldenrod1","dodgerblue3")
color_possibilities_RE <- c("Control","UC","CD")
names(myColors_RE) <- color_possibilities_RE
fillscale_RE_diagnosis2 <- scale_fill_manual(name="Diagnosis", values=myColors_RE, drop=T, limits=force)
colscale_RE_diagnosis2 <- scale_color_manual(name="Diagnosis", values=myColors_RE, drop=T)

myColours <- c("grey30","grey70")
colour_possibilities <- c("severe","other")
names(myColours) <- colour_possibilities
fillscale_2_severity <- scale_fill_manual(name="Disease Severity", labels=c("severe"="severe","other"="mild/moderate"), values=myColours, drop=T, limits=force, na.value="white")
colscale_2_severity <- scale_color_manual(name="Disease Severity", values=myColours, drop=T)


######################################
# HEAT SCREE PLOTS: MHCI PROBES ONLY #
######################################

# ALL SAMPLES:

ti.pca <- prcomp(t(ti.mhc))
Loadings <- as.data.frame(ti.pca$x)
vars <- ti.pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- ti.samples[c("sex","TI.numeric","Biologics.numeric","Surgery.numeric","AZA.numeric","Perianal.numeric","TI.MHCI.Group")]
meta_categorical$TI.MHCI.Group <- ifelse(meta_categorical$TI.MHCI.Group=="Low",1,0)
meta_continuous <- ti.samples[c("age","Severity")]
colnames(meta_categorical) <- c("Sex","Inflammation","Biologics","Surgery","AZA","Perianal Disease","MHCI Group")
colnames(meta_continuous) <- c("Age","Disease Severity")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_TI_MHCIProbes_OutcomeMeasures_HeatScree_31Oct22.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

sc.pca <- prcomp(t(sc.mhc))
Loadings <- as.data.frame(sc.pca$x)
vars <- ti.pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- sc.samples[c("sex","SC.numeric","Biologics.numeric","Surgery.numeric","AZA.numeric","Perianal.numeric","SC.MHCI.Group")]
meta_categorical$SC.MHCI.Group <- ifelse(meta_categorical$SC.MHCI.Group=="Low",1,0)
meta_continuous <- sc.samples[c("age","Severity")]
colnames(meta_categorical) <- c("Sex","Inflammation","Biologics","Surgery","AZA","Perianal Disease","MHCI Group")
colnames(meta_continuous) <- c("Age","Disease Severity")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_SC_MHCIProbes_OutcomeMeasures_HeatScree_31Oct22.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()


# CD SAMPLES:

cd.samples <- subset(ti.samples, diagnosis=="CD")
cd.mhc <- ti.mhc[,cd.samples$array.id]
ti.pca <- prcomp(t(cd.mhc))
Loadings <- as.data.frame(ti.pca$x)
vars <- ti.pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- cd.samples[c("sex","TI.numeric","Biologics.numeric","Surgery.numeric","AZA.numeric","Perianal.numeric","TI.MHCI.Group")]
meta_categorical$TI.MHCI.Group <- ifelse(meta_categorical$TI.MHCI.Group=="Low",1,0)
meta_continuous <- cd.samples[c("age","Severity")]
colnames(meta_categorical) <- c("Sex","Inflammation","Biologics","Surgery","AZA","Perianal Disease","MHCI Group")
colnames(meta_continuous) <- c("Age","Disease Severity")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_TI-CD_MHCIProbes_OutcomeMeasures_HeatScree_31Oct22.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

cd.sc.samples <- subset(sc.samples, diagnosis=="CD")
cd.sc.mhc <- sc.mhc[,cd.sc.samples$array.id]
sc.pca <- prcomp(t(cd.sc.mhc))
Loadings <- as.data.frame(sc.pca$x)
vars <- sc.pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- cd.sc.samples[c("sex","SC.numeric","Biologics.numeric","Surgery.numeric","AZA.numeric","Perianal.numeric","SC.MHCI.Group")]
meta_categorical$SC.MHCI.Group <- ifelse(meta_categorical$SC.MHCI.Group=="Low",1,0)
meta_continuous <- cd.sc.samples[c("age","Severity")]
colnames(meta_categorical) <- c("Sex","Inflammation","Biologics","Surgery","AZA","Perianal Disease","MHCI Group")
colnames(meta_continuous) <- c("Age","Disease Severity")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_SC-CD_MHCIProbes_OutcomeMeasures_HeatScree_31Oct22.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()


# IBD SAMPLES:

ibd.samples <- subset(ti.samples, diagnosis!="Control")
ibd.mhc <- ti.mhc[,ibd.samples$array.id]
ti.pca <- prcomp(t(ibd.mhc))
Loadings <- as.data.frame(ti.pca$x)
vars <- ti.pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- ibd.samples[c("sex","TI.numeric","Biologics.numeric","Surgery.numeric","AZA.numeric","Perianal.numeric","TI.MHCI.Group")]
meta_categorical$TI.MHCI.Group <- ifelse(meta_categorical$TI.MHCI.Group=="Low",1,0)
meta_continuous <- ibd.samples[c("age","Severity")]
colnames(meta_categorical) <- c("Sex","Inflammation","Biologics","Surgery","AZA","Perianal Disease","MHCI Group")
colnames(meta_continuous) <- c("Age","Disease Severity")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_TI-IBD_MHCIProbes_OutcomeMeasures_HeatScree_31Oct22.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()

ibd.sc.samples <- subset(sc.samples, diagnosis!="Control")
ibd.sc.mhc <- sc.mhc[,ibd.sc.samples$array.id]
sc.pca <- prcomp(t(ibd.sc.mhc))
Loadings <- as.data.frame(sc.pca$x)
vars <- sc.pca$sdev^2
Importance <- vars/sum(vars)

meta_categorical <- ibd.sc.samples[c("sex","SC.numeric","Biologics.numeric","Surgery.numeric","AZA.numeric","Perianal.numeric","SC.MHCI.Group")]
meta_categorical$SC.MHCI.Group <- ifelse(meta_categorical$SC.MHCI.Group=="Low",1,0)
meta_continuous <- ibd.sc.samples[c("age","Severity")]
colnames(meta_categorical) <- c("Sex","Inflammation","Biologics","Surgery","AZA","Perianal Disease","MHCI Group")
colnames(meta_continuous) <- c("Age","Disease Severity")
ord <- 1:length(c(colnames(meta_categorical),colnames(meta_continuous)))

pdf("HEAT_SCREE_PLOTS/PaediatricOrganoid_SC-IBD_MHCIProbes_OutcomeMeasures_HeatScree_31Oct22.pdf", width=10, height=6)
heat.scree(Loadings, Importance, 2.5, 2.7)
dev.off()


#####################################################
# AVERAGE MHCI METHYLATION SCORE DISTRIBUTION PLOTS #
#####################################################

# 3 CLUSTERS:
# Groups: Low/Intermediate/High

ti.distrib <- as.data.frame(table(clust$TI.MHCI.Group.2,clust$diagnosis))
colnames(ti.distrib) <- c("Cluster","Diagnosis","Number")
sc.distrib <- as.data.frame(table(clust$SC.MHCI.Group.2,clust$diagnosis))
colnames(sc.distrib) <- c("Cluster","Diagnosis","Number")

# Percentage stacked barplots of sample distribution by diagnosis:

p1 <- ggplot(ti.distrib, aes(fill=Cluster, y=Number, x=Diagnosis))+geom_bar(position="fill", stat="identity")

pdf("MHCI_SCORE_DISTRIBUTION/PaediatricOrganoid_TI_SampleDistributionByCluster_3Clusters_31Oct22.pdf")
p1+theme_bw()+fillscale_3_cluster+ylab("Percentage of Samples (TI)")+xlab("")+scale_x_discrete(limits=c("CD","UC","Control"))
dev.off()

p2 <- ggplot(sc.distrib, aes(fill=Cluster, y=Number, x=Diagnosis))+geom_bar(position="fill", stat="identity")

pdf("MHCI_SCORE_DISTRIBUTION/PaediatricOrganoid_SC_SampleDistributionByCluster_3Clusters_31Oct22.pdf")
p2+theme_bw()+fillscale_3_cluster+ylab("Percentage of Samples (SC)")+xlab("")+scale_x_discrete(limits=c("CD","UC","Control"))
dev.off()

# Percentage stacked barplots of sample distribution by MHCI group:

p1 <- ggplot(ti.distrib, aes(fill=Diagnosis, y=Number, x=Cluster))+geom_bar(position="fill", stat="identity")

pdf("MHCI_SCORE_DISTRIBUTION/PaediatricOrganoid_TI_SampleDistributionByDiagnosis_3Clusters_31Oct22.pdf")
p1+theme_bw()+fillscale_RE_diagnosis2+ylab("Percentage of Samples (TI)")+xlab("MHCI Average Methylation Group")+scale_x_discrete(limits=c("Low","Intermediate","High"))
dev.off()

p2 <- ggplot(sc.distrib, aes(fill=Diagnosis, y=Number, x=Cluster))+geom_bar(position="fill", stat="identity")

pdf("MHCI_SCORE_DISTRIBUTION/PaediatricOrganoid_SC_SampleDistributionByDiagnosis_3Clusters_31Oct22.pdf")
p2+theme_bw()+fillscale_RE_diagnosis2+ylab("Percentage of Samples (SC)")+xlab("MHCI Average Methylation Group")+scale_x_discrete(limits=c("Low","Intermediate","High"))
dev.off()


# 2 CLUSTERS:
# Groups: Low/Other

ti.distrib <- as.data.frame(table(clust$TI.MHCI.Group,clust$diagnosis))
colnames(ti.distrib) <- c("Cluster","Diagnosis","Number")

sc.distrib <- as.data.frame(table(clust$SC.MHCI.Group,clust$diagnosis))
colnames(sc.distrib) <- c("Cluster","Diagnosis","Number")

# Percentage stacked barplots of sample distribution by diagnosis:

p1 <- ggplot(ti.distrib, aes(fill=Cluster, y=Number, x=Diagnosis))+geom_bar(position="fill", stat="identity")

pdf("MHCI_SCORE_DISTRIBUTION/PaediatricOrganoid_TI_SampleDistributionByCluster_2Clusters_31Oct22.pdf")
p1+theme_bw()+fillscale_2_cluster+ylab("Percentage of Samples (TI)")+xlab("")+scale_x_discrete(limits=c("CD","UC","Control"))
dev.off()

p2 <- ggplot(sc.distrib, aes(fill=Cluster, y=Number, x=Diagnosis))+geom_bar(position="fill", stat="identity")

pdf("MHCI_SCORE_DISTRIBUTION/PaediatricOrganoid_SC_SampleDistributionByCluster_2Clusters_31Oct22.pdf")
p2+theme_bw()+fillscale_2_cluster+ylab("Percentage of Samples (SC)")+xlab("")+scale_x_discrete(limits=c("CD","UC","Control"))
dev.off()

# Percentage stacked barplots of sample distribution by MHCI group:

p1 <- ggplot(ti.distrib, aes(fill=Diagnosis, y=Number, x=Cluster))+geom_bar(position="fill", stat="identity")

pdf("MHCI_SCORE_DISTRIBUTION/PaediatricOrganoid_TI_SampleDistributionByDiagnosis_2Clusters_31Oct22.pdf")
p1+theme_bw()+fillscale_RE_diagnosis2+ylab("Percentage of Samples (TI)")+xlab("MHCI Average Methylation Group")+scale_x_discrete(limits=c("Low","Other"))
dev.off()

p2 <- ggplot(sc.distrib, aes(fill=Diagnosis, y=Number, x=Cluster))+geom_bar(position="fill", stat="identity")

pdf("MHCI_SCORE_DISTRIBUTION/PaediatricOrganoid_SC_SampleDistributionByDiagnosis_2Clusters_31Oct22.pdf")
p2+theme_bw()+fillscale_RE_diagnosis2+ylab("Percentage of Samples (SC)")+xlab("MHCI Average Methylation Group")+scale_x_discrete(limits=c("Low","Other"))
dev.off()


##################################################
# AVERAGE METHYLATION BY MHCI GROUP (2 CLUSTERS) #
##################################################

# Groups: Low/Other
# NB, there are no Control/UC/IBDU SC samples in the "Low" group, so I need to add (non-significant) false data points for these in order to plot using the group_by function

ti.avg <- ti.samples[,c("array.id","diagnosis","TI.MHCI.Group","Avg.MHCI.CpG")]
ti.avg <- melt(ti.avg, c("array.id","diagnosis","TI.MHCI.Group"), na.rm=T)
colnames(ti.avg) <- c("Sample","Diagnosis","Group","CpG","Beta")

sc.avg <- sc.samples[,c("array.id","diagnosis","SC.MHCI.Group","Avg.MHCI.CpG")]
sc.avg <- melt(sc.avg, c("array.id","diagnosis","SC.MHCI.Group"), na.rm=T)
colnames(sc.avg) <- c("Sample","Diagnosis","Group","CpG","Beta")

table(ti.avg$Group,ti.avg$Diagnosis)
table(sc.avg$Group,sc.avg$Diagnosis)

# Plot average MHCI methylation score by diagnosis:

p1 <- ggplot(ti.avg, aes(x=Group, y=Beta, group=Group))+geom_boxplot(aes(fill=Group), width=0.7)
s1 <- ti.avg %>% dplyr::group_by(Diagnosis) %>% wilcox_test(Beta~Group) %>% adjust_pvalue(method="BH") %>% add_significance()
s1 <- s1 %>% add_y_position()

s1

s1$p.sci <- ifelse(s1$p.adj<=0.05, (formatC(s1$p.adj, digits=2, format="e")), "ns")

pdf("MHCI_SCORE_DISTRIBUTION/PaediatricOrganoid_TI_AvgMHCI_ByDiagnosis_31Oct22.pdf", height=5, width=11)
p1+facet_grid(cols=vars(Diagnosis))+theme_bw()+scale_x_discrete(limits=c("Low","Other"))+xlab("MHCI Average Methylation Group")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Group),color="black",shape=21,size=2,alpha=0.9,width=0.2)+fillscale_2_cluster+stat_pvalue_manual(s1, hide.ns=T, label="p.sci", tip.length=0.01)
dev.off()

p2 <- ggplot(sc.avg, aes(x=Group, y=Beta, group=Group))+geom_boxplot(aes(fill=Group), width=0.7)
false.cont <- mean((subset(sc.avg, Diagnosis=="Control"))$Beta)
false.uc <- mean((subset(sc.avg, Diagnosis=="UC"))$Beta)
false <- data.frame(Sample=c("false1","false2"),Diagnosis=c("Control","UC"),Group=c("Low","Low"),CpG=c("Avg.MHCI.CpG","Avg.MHCI.CpG"),Beta=c(false.cont,false.uc))
false <- rbind(false,sc.avg)
s2 <- false %>% dplyr::group_by(Diagnosis) %>% wilcox_test(Beta~Group) %>% adjust_pvalue(method="BH") %>% add_significance()
s2 <- s2 %>% add_y_position()

s2

s2$p.sci <- ifelse(s2$p.adj<=0.05, (formatC(s2$p.adj, digits=2, format="e")), "ns")

pdf("MHCI_SCORE_DISTRIBUTION/PaediatricOrganoid_SC_AvgMHCI_ByDiagnosis_31Oct22.pdf", height=5, width=11)
p2+facet_grid(cols=vars(Diagnosis))+theme_bw()+scale_x_discrete(limits=c("Low","Other"))+xlab("MHCI Average Methylation Group")+ylab("Average MHCI Methylation (SC)")+geom_jitter(aes(fill=Group),color="black",shape=21,size=2,alpha=0.9,width=0.2)+fillscale_2_cluster+stat_pvalue_manual(s2, hide.ns=T, label="p.sci", tip.length=0.01)
dev.off()


#########################################################################
# CLINICAL DATA CORRELATIONS (2 CLUSTERS): CD ONLY - FISHER GROUP TESTS #
#########################################################################

# TI/SC CD samples only
# Groups: Low/Other

library(viridis)

# TI SAMPLES:

cd.samples <- subset(ti.samples, diagnosis=="CD")
cd.samples.2 <- subset(cd.samples, Disease.Severity!="moderate")

# TI Inflammation:

table(cd.samples$TI.MHCI.Group,cd.samples$TI.numeric)

f1 <- fisher.test(table(cd.samples$TI.MHCI.Group,cd.samples$TI.numeric))
pval <- formatC(f1$p.value,digits=2,format="e")
ti.fisher <- data.frame(Feature="Inflammation",Low.v.Other.Fisher=pval)

f1

# Biologics:

table(cd.samples$TI.MHCI.Group,cd.samples$Biologics.numeric)

f2 <- fisher.test(table(cd.samples$TI.MHCI.Group,cd.samples$Biologics.numeric))
pval <- formatC(f2$p.value,digits=2,format="e")
update <- data.frame(Feature="Biologics",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f2

# Surgery:

table(cd.samples$TI.MHCI.Group,cd.samples$Surgery.numeric)

f3 <- fisher.test(table(cd.samples$TI.MHCI.Group,cd.samples$Surgery.numeric))
pval <- formatC(f3$p.value,digits=2,format="e")
update <- data.frame(Feature="Surgery",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f3

# AZA treatment:

table(cd.samples$TI.MHCI.Group,cd.samples$AZA.numeric)

f4 <- fisher.test(table(cd.samples$TI.MHCI.Group,cd.samples$AZA.numeric))
pval <- formatC(f4$p.value,digits=2,format="e")
update <- data.frame(Feature="AZA Treatment",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f4

# Perianal disease:

table(cd.samples$TI.MHCI.Group,cd.samples$Perianal.numeric)

f5 <- fisher.test(table(cd.samples$TI.MHCI.Group,cd.samples$Perianal.numeric))
pval <- formatC(f5$p.value,digits=2,format="e")
update <- data.frame(Feature="Perianal Disease",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f5

# Disease severity: severe vs mild/moderate

table(cd.samples$TI.MHCI.Group,cd.samples$severe.v.other)

f6 <- fisher.test(table(cd.samples$TI.MHCI.Group,cd.samples$severe.v.other))
pval <- formatC(f6$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe v Mild/Moderate",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f6

# Disease severity: severe/moderate vs mild

table(cd.samples$TI.MHCI.Group,cd.samples$mild.v.other)

f7 <- fisher.test(table(cd.samples$TI.MHCI.Group,cd.samples$mild.v.other))
pval <- formatC(f7$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe/Moderate v Mild",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f7

# Disease severity: severe vs mild

table(cd.samples.2$TI.MHCI.Group,cd.samples.2$mild.v.severe)

f8 <- fisher.test(table(cd.samples.2$TI.MHCI.Group,cd.samples.2$Disease.Severity))
pval <- formatC(f8$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe v Mild (subset)",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f8


# SC SAMPLES:

cd.sc.samples <- subset(sc.samples, diagnosis=="CD")
cd.sc.samples.2 <- subset(cd.sc.samples, Disease.Severity!="moderate")

# SC Inflammation:

table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$SC.numeric)

f1 <- fisher.test(table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$SC.numeric))
pval <- formatC(f1$p.value,digits=2,format="e")
sc.fisher <- data.frame(Feature="Inflammation",Low.v.Other.Fisher=pval)

f1

# Biologics:

table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$Biologics.numeric)

f2 <- fisher.test(table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$Biologics.numeric))
pval <- formatC(f2$p.value,digits=2,format="e")
update <- data.frame(Feature="Biologics",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f2

# Surgery:

table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$Surgery.numeric)

f3 <- fisher.test(table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$Surgery.numeric))
pval <- formatC(f3$p.value,digits=2,format="e")
update <- data.frame(Feature="Surgery",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f3

# AZA treatment:

table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$AZA.numeric)

f4 <- fisher.test(table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$AZA.numeric))
pval <- formatC(f4$p.value,digits=2,format="e")
update <- data.frame(Feature="AZA Treatment",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f4

# Perianal disease:

table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$Perianal.numeric)

f5 <- fisher.test(table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$Perianal.numeric))
pval <- formatC(f5$p.value,digits=2,format="e")
update <- data.frame(Feature="Perianal Disease",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f5

# Disease severity: severe vs mild/moderate

table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$severe.v.other)

f6 <- fisher.test(table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$severe.v.other))
pval <- formatC(f6$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe v Mild/Moderate",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f6

# Disease severity: severe/moderate vs mild

table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$mild.v.other)

f7 <- fisher.test(table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$mild.v.other))
pval <- formatC(f7$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe/Moderate v Mild",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f7

# Disease severity: severe vs mild

table(cd.sc.samples.2$SC.MHCI.Group,cd.sc.samples.2$mild.v.severe)

f8 <- fisher.test(table(cd.sc.samples.2$SC.MHCI.Group,cd.sc.samples.2$Disease.Severity))
pval <- formatC(f8$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe v Mild (subset)",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f8


#########################################################################
# CLINICAL DATA CORRELATIONS (2 CLUSTERS): CD ONLY - WILCOX SCORE TESTS #
#########################################################################

# CLINICAL DATA CORRELATIONS (Average MHCI Score):
# TI/SC CD samples only

# TI SAMPLES:
# TI Inflammation:

cd.samples$TI.Inflammation.2 <- ifelse(cd.samples$TI.numeric==1,"Y","N")

w1 <- wilcox_test(cd.samples, Avg.MHCI.CpG~TI.numeric, p.adjust.method="BH")

w1

w1 <- w1 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-CD_AvgMHCI_Inflammation_31Oct22.pdf")
ggplot(cd.samples, aes(x=TI.Inflammation.2, y=Avg.MHCI.CpG, group=TI.Inflammation.2))+geom_boxplot(aes(fill=TI.Inflammation.2), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("TI Inflammation")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=TI.Inflammation.2), color="black",shape=21,size=2,alpha=0.9)+fillscale_inflammation2+stat_pvalue_manual(w1, label="p", tip.length=0.01)
dev.off()

ti.wilcox <- data.frame(Feature="Inflammation",AvgMHC.Wilcox=w1$p)

# Biologics:

w2 <- wilcox_test(cd.samples, Avg.MHCI.CpG~Biologics.numeric, p.adjust.method="BH")

w2

w2 <- w2 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-CD_AvgMHCI_Biologics_31Oct22.pdf")
ggplot(cd.samples, aes(x=Biologics, y=Avg.MHCI.CpG, group=Biologics))+geom_boxplot(aes(fill=Biologics), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Biologics")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Biologics), color="black",shape=21,size=2,alpha=0.9)+fillscale_biologics+stat_pvalue_manual(w2, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Biologics",AvgMHC.Wilcox=w2$p)
ti.wilcox <- rbind(ti.wilcox,update)

# Surgery:

w3 <- wilcox_test(cd.samples, Avg.MHCI.CpG~Surgery.numeric, p.adjust.method="BH")

w3

w3 <- w3 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-CD_AvgMHCI_Surgery_31Oct22.pdf")
ggplot(cd.samples, aes(x=Surgery, y=Avg.MHCI.CpG, group=Surgery))+geom_boxplot(aes(fill=Surgery), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Surgery")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Surgery), color="black",shape=21,size=2,alpha=0.9)+fillscale_surgery+stat_pvalue_manual(w3, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Surgery",AvgMHC.Wilcox=w3$p)
ti.wilcox <- rbind(ti.wilcox,update)


# AZA treatment:

w4 <- wilcox_test(cd.samples, Avg.MHCI.CpG~AZA.numeric, p.adjust.method="BH")

w4

w4 <- w4 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-CD_AvgMHCI_AZATreatment_31Oct22.pdf")
ggplot(cd.samples, aes(x=AZA, y=Avg.MHCI.CpG, group=AZA))+geom_boxplot(aes(fill=AZA), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("AZA")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=AZA), color="black",shape=21,size=2,alpha=0.9)+fillscale_immunosuppressor+stat_pvalue_manual(w4, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="AZA Treatment",AvgMHC.Wilcox=w4$p)
ti.wilcox <- rbind(ti.wilcox,update)

# Perianal disease:

w5 <- wilcox_test(cd.samples, Avg.MHCI.CpG~Perianal.numeric, p.adjust.method="BH")

w5

w5 <- w5 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-CD_AvgMHCI_PerianalDisease_31Oct22.pdf")
ggplot(cd.samples, aes(x=Perianal.Disease, y=Avg.MHCI.CpG, group=Perianal.Disease))+geom_boxplot(aes(fill=Perianal.Disease), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Perianal Disease")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Perianal.Disease), color="black",shape=21,size=2,alpha=0.9)+fillscale_perianal+stat_pvalue_manual(w5, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Perianal Disease",AvgMHC.Wilcox=w5$p)
ti.wilcox <- rbind(ti.wilcox,update)

# Disease severity: severe vs mild/moderate

w6 <- wilcox_test(cd.samples, Avg.MHCI.CpG~severe.v.other, p.adjust.method="BH")

w6

w6 <- w6 %>% add_y_position()
cd.samples$severe.v.other <- as.character(cd.samples$severe.v.other)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-CD_AvgMHCI_Severe-Other_31Oct22.pdf")
ggplot(cd.samples, aes(x=severe.v.other, y=Avg.MHCI.CpG, group=severe.v.other))+geom_boxplot(aes(fill=severe.v.other), width=0.7)+theme_bw()+scale_x_discrete(limits=c("1","0"), labels=c("1"="Severe","0"="Mild/Moderate"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=severe.v.other), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity2+stat_pvalue_manual(w6, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe v Mild/Moderate",AvgMHC.Wilcox=w6$p)
ti.wilcox <- rbind(ti.wilcox,update)

# Disease severity: severe/moderate vs mild

w7 <- wilcox_test(cd.samples, Avg.MHCI.CpG~mild.v.other, p.adjust.method="BH")

w7

w7 <- w7 %>% add_y_position()
cd.samples$mild.v.other <- as.character(cd.samples$mild.v.other)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-CD_AvgMHCI_Mild-Other_31Oct22.pdf")
ggplot(cd.samples, aes(x=mild.v.other, y=Avg.MHCI.CpG, group=mild.v.other))+geom_boxplot(aes(fill=mild.v.other), width=0.7)+theme_bw()+scale_x_discrete(limits=c("0","1"), labels=c("0"="Moderate/Severe","1"="Mild"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=mild.v.other), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity3+stat_pvalue_manual(w7, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe/Moderate v Mild",AvgMHC.Wilcox=w7$p)
ti.wilcox <- rbind(ti.wilcox,update)

# Disease severity: severe vs mild

w8 <- wilcox_test(cd.samples.2, Avg.MHCI.CpG~mild.v.severe, p.adjust.method="BH")

w8

w8 <- w8 %>% add_y_position()
cd.samples.2$mild.v.severe <- as.character(cd.samples.2$mild.v.severe)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-CD_AvgMHCI_Mild-Severe_31Oct22.pdf")
ggplot(cd.samples.2, aes(x=mild.v.severe, y=Avg.MHCI.CpG, group=mild.v.severe))+geom_boxplot(aes(fill=mild.v.severe), width=0.7)+theme_bw()+scale_x_discrete(limits=c("1","0"), labels=c("1"="Severe","0"="Mild"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=mild.v.severe), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity2+stat_pvalue_manual(w8, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe v Mild (subset)",AvgMHC.Wilcox=w8$p)
ti.wilcox <- rbind(ti.wilcox,update)

# Write all TI fisher's exact & wilcox test results to file:

ti.stats <- join(ti.wilcox, ti.fisher, type="left", match="all")
write.table(ti.stats, "MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-CD_MHCI_Correlations_31Oct22.txt", row.names=F, sep="\t", quote=F)


# SC SAMPLES:
# SC Inflammation:

cd.sc.samples$SC.Inflammation.2 <- ifelse(cd.sc.samples$SC.numeric==1,"Y","N")

w1 <- wilcox_test(cd.sc.samples, Avg.MHCI.CpG~SC.numeric, p.adjust.method="BH")

w1

w1 <- w1 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-CD_AvgMHCI_Inflammation_31Oct22.pdf")
ggplot(cd.sc.samples, aes(x=SC.Inflammation.2, y=Avg.MHCI.CpG, group=SC.Inflammation.2))+geom_boxplot(aes(fill=SC.Inflammation.2), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("SC Inflammation")+ylab("Average MHCI Methylation (SC)")+geom_jitter(aes(fill=SC.Inflammation.2), color="black",shape=21,size=2,alpha=0.9)+fillscale_inflammation2+stat_pvalue_manual(w1, label="p", tip.length=0.01)
dev.off()

sc.wilcox <- data.frame(Feature="Inflammation",AvgMHC.Wilcox=w1$p)

# Biologics:

w2 <- wilcox_test(cd.sc.samples, Avg.MHCI.CpG~Biologics.numeric, p.adjust.method="BH")

w2

w2 <- w2 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-CD_AvgMHCI_Biologics_31Oct22.pdf")
ggplot(cd.sc.samples, aes(x=Biologics, y=Avg.MHCI.CpG, group=Biologics))+geom_boxplot(aes(fill=Biologics), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Biologics")+ylab("Average MHCI Methylation (SC)")+geom_jitter(aes(fill=Biologics), color="black",shape=21,size=2,alpha=0.9)+fillscale_biologics+stat_pvalue_manual(w2, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Biologics",AvgMHC.Wilcox=w2$p)
sc.wilcox <- rbind(sc.wilcox,update)

# Surgery:

w3 <- wilcox_test(cd.sc.samples, Avg.MHCI.CpG~Surgery.numeric, p.adjust.method="BH")

w3

w3 <- w3 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-CD_AvgMHCI_Surgery_31Oct22.pdf")
ggplot(cd.sc.samples, aes(x=Surgery, y=Avg.MHCI.CpG, group=Surgery))+geom_boxplot(aes(fill=Surgery), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Surgery")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Surgery), color="black",shape=21,size=2,alpha=0.9)+fillscale_surgery+stat_pvalue_manual(w3, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Surgery",AvgMHC.Wilcox=w3$p)
sc.wilcox <- rbind(sc.wilcox,update)


# AZA treatment:

w4 <- wilcox_test(cd.sc.samples, Avg.MHCI.CpG~AZA.numeric, p.adjust.method="BH")

w4

w4 <- w4 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-CD_AvgMHCI_AZATreatment_31Oct22.pdf")
ggplot(cd.sc.samples, aes(x=AZA, y=Avg.MHCI.CpG, group=AZA))+geom_boxplot(aes(fill=AZA), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("AZA")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=AZA), color="black",shape=21,size=2,alpha=0.9)+fillscale_immunosuppressor+stat_pvalue_manual(w4, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="AZA Treatment",AvgMHC.Wilcox=w4$p)
sc.wilcox <- rbind(sc.wilcox,update)

# Perianal disease:

w5 <- wilcox_test(cd.sc.samples, Avg.MHCI.CpG~Perianal.numeric, p.adjust.method="BH")

w5

w5 <- w5 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-CD_AvgMHCI_PerianalDisease_31Oct22.pdf")
ggplot(cd.sc.samples, aes(x=Perianal.Disease, y=Avg.MHCI.CpG, group=Perianal.Disease))+geom_boxplot(aes(fill=Perianal.Disease), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Perianal Disease")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Perianal.Disease), color="black",shape=21,size=2,alpha=0.9)+fillscale_perianal+stat_pvalue_manual(w5, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Perianal Disease",AvgMHC.Wilcox=w5$p)
sc.wilcox <- rbind(sc.wilcox,update)

# Disease severity: severe vs mild/moderate

w6 <- wilcox_test(cd.sc.samples, Avg.MHCI.CpG~severe.v.other, p.adjust.method="BH")

w6

w6 <- w6 %>% add_y_position()
cd.sc.samples$severe.v.other <- as.character(cd.sc.samples$severe.v.other)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-CD_AvgMHCI_Severe-Other_31Oct22.pdf")
ggplot(cd.sc.samples, aes(x=severe.v.other, y=Avg.MHCI.CpG, group=severe.v.other))+geom_boxplot(aes(fill=severe.v.other), width=0.7)+theme_bw()+scale_x_discrete(limits=c("1","0"), labels=c("1"="Severe","0"="Mild/Moderate"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=severe.v.other), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity2+stat_pvalue_manual(w6, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe v Mild/Moderate",AvgMHC.Wilcox=w6$p)
sc.wilcox <- rbind(sc.wilcox,update)

# Disease severity: severe/moderate vs mild

w7 <- wilcox_test(cd.sc.samples, Avg.MHCI.CpG~mild.v.other, p.adjust.method="BH")

w7

w7 <- w7 %>% add_y_position()
cd.sc.samples$mild.v.other <- as.character(cd.sc.samples$mild.v.other)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-CD_AvgMHCI_Mild-Other_31Oct22.pdf")
ggplot(cd.sc.samples, aes(x=mild.v.other, y=Avg.MHCI.CpG, group=mild.v.other))+geom_boxplot(aes(fill=mild.v.other), width=0.7)+theme_bw()+scale_x_discrete(limits=c("0","1"), labels=c("0"="Moderate/Severe","1"="Mild"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=mild.v.other), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity3+stat_pvalue_manual(w7, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe/Moderate v Mild",AvgMHC.Wilcox=w7$p)
sc.wilcox <- rbind(sc.wilcox,update)

# Disease severity: severe vs mild

w8 <- wilcox_test(cd.sc.samples.2, Avg.MHCI.CpG~mild.v.severe, p.adjust.method="BH")

w8

w8 <- w8 %>% add_y_position()
cd.sc.samples.2$mild.v.severe <- as.character(cd.sc.samples.2$mild.v.severe)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-CD_AvgMHCI_Mild-Severe_31Oct22.pdf")
ggplot(cd.sc.samples.2, aes(x=mild.v.severe, y=Avg.MHCI.CpG, group=mild.v.severe))+geom_boxplot(aes(fill=mild.v.severe), width=0.7)+theme_bw()+scale_x_discrete(limits=c("1","0"), labels=c("1"="Severe","0"="Mild"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=mild.v.severe), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity2+stat_pvalue_manual(w8, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe v Mild (subset)",AvgMHC.Wilcox=w8$p)
sc.wilcox <- rbind(sc.wilcox,update)

# Write all SC fisher's exact & wilcox test results to file:

sc.stats <- join(sc.wilcox, sc.fisher, type="left", match="all")
write.table(sc.stats, "MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-CD_MHCI_Correlations_31Oct22.txt", row.names=F, sep="\t", quote=F)


##########################################################################
# CLINICAL DATA CORRELATIONS (2 CLUSTERS): IBD ONLY - FISHER GROUP TESTS #
##########################################################################

# TI/SC IBD samples only (ie everything except the controls)
# Groups: Low/Other

# TI SAMPLES:

ibd.samples <- subset(ti.samples, diagnosis!="Control")
ibd.samples$severe.v.other <- ifelse(ibd.samples$Disease.Severity=="severe",1,0)
ibd.samples$mild.v.other <- ifelse(ibd.samples$Disease.Severity=="mild",0,1)
ibd.samples.2 <- subset(ibd.samples, Disease.Severity!="moderate")
ibd.samples.2$mild.v.severe <- ifelse(ibd.samples.2$Disease.Severity=="severe",1,0)

# TI Inflammation:

table(ibd.samples$TI.MHCI.Group,ibd.samples$TI.numeric)

f1 <- fisher.test(table(ibd.samples$TI.MHCI.Group,ibd.samples$TI.numeric))
pval <- formatC(f1$p.value,digits=2,format="e")
ti.fisher <- data.frame(Feature="Inflammation",Low.v.Other.Fisher=pval)

f1

# Biologics:

table(ibd.samples$TI.MHCI.Group,ibd.samples$Biologics.numeric)

f2 <- fisher.test(table(ibd.samples$TI.MHCI.Group,ibd.samples$Biologics.numeric))
pval <- formatC(f2$p.value,digits=2,format="e")
update <- data.frame(Feature="Biologics",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f2

# Surgery:

table(ibd.samples$TI.MHCI.Group,ibd.samples$Surgery.numeric)

f3 <- fisher.test(table(ibd.samples$TI.MHCI.Group,ibd.samples$Surgery.numeric))
pval <- formatC(f3$p.value,digits=2,format="e")
update <- data.frame(Feature="Surgery",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f3

# AZA treatment:

table(ibd.samples$TI.MHCI.Group,ibd.samples$AZA.numeric)

f4 <- fisher.test(table(ibd.samples$TI.MHCI.Group,ibd.samples$AZA.numeric))
pval <- formatC(f4$p.value,digits=2,format="e")
update <- data.frame(Feature="AZA Treatment",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f4

# Perianal disease:

table(ibd.samples$TI.MHCI.Group,ibd.samples$Perianal.numeric)

f5 <- fisher.test(table(ibd.samples$TI.MHCI.Group,ibd.samples$Perianal.numeric))
pval <- formatC(f5$p.value,digits=2,format="e")
update <- data.frame(Feature="Perianal Disease",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f5

# Disease severity: severe vs mild/moderate

table(ibd.samples$TI.MHCI.Group,ibd.samples$severe.v.other)

f6 <- fisher.test(table(ibd.samples$TI.MHCI.Group,ibd.samples$severe.v.other))
pval <- formatC(f6$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe v Mild/Moderate",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f6

# Disease severity: severe/moderate vs mild

table(ibd.samples$TI.MHCI.Group,ibd.samples$mild.v.other)

f7 <- fisher.test(table(ibd.samples$TI.MHCI.Group,ibd.samples$mild.v.other))
pval <- formatC(f7$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe/Moderate v Mild",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f7

# Disease severity: severe vs mild

table(ibd.samples.2$TI.MHCI.Group,ibd.samples.2$mild.v.severe)

f8 <- fisher.test(table(ibd.samples.2$TI.MHCI.Group,ibd.samples.2$mild.v.severe))
pval <- formatC(f8$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe v Mild (subset)",Low.v.Other.Fisher=pval)
ti.fisher <- rbind(ti.fisher,update)

f8


# SC SAMPLES:

ibd.sc.samples <- subset(sc.samples, diagnosis!="Control")
ibd.sc.samples$severe.v.other <- ifelse(ibd.sc.samples$Disease.Severity=="severe",1,0)
ibd.sc.samples$mild.v.other <- ifelse(ibd.sc.samples$Disease.Severity=="mild",0,1)
ibd.sc.samples.2 <- subset(ibd.sc.samples, Disease.Severity!="moderate")
ibd.sc.samples.2$mild.v.severe <- ifelse(ibd.sc.samples.2$Disease.Severity=="severe",1,0)

# SC Inflammation:

table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$SC.numeric)

f1 <- fisher.test(table(cd.sc.samples$SC.MHCI.Group,cd.sc.samples$SC.numeric))
pval <- formatC(f1$p.value,digits=2,format="e")
sc.fisher <- data.frame(Feature="Inflammation",Low.v.Other.Fisher=pval)

f1

# Biologics:

table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$Biologics.numeric)

f2 <- fisher.test(table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$Biologics.numeric))
pval <- formatC(f2$p.value,digits=2,format="e")
update <- data.frame(Feature="Biologics",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f2

# Surgery:

table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$Surgery.numeric)

f3 <- fisher.test(table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$Surgery.numeric))
pval <- formatC(f3$p.value,digits=2,format="e")
update <- data.frame(Feature="Surgery",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f3

# AZA treatment:

table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$AZA.numeric)

f4 <- fisher.test(table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$AZA.numeric))
pval <- formatC(f4$p.value,digits=2,format="e")
update <- data.frame(Feature="AZA Treatment",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f4

# Perianal disease:

table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$Perianal.numeric)

f5 <- fisher.test(table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$Perianal.numeric))
pval <- formatC(f5$p.value,digits=2,format="e")
update <- data.frame(Feature="Perianal Disease",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f5

# Disease severity: severe vs mild/moderate

table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$severe.v.other)

f6 <- fisher.test(table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$severe.v.other))
pval <- formatC(f6$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe v Mild/Moderate",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f6

# Disease severity: severe/moderate vs mild

table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$mild.v.other)

f7 <- fisher.test(table(ibd.sc.samples$SC.MHCI.Group,ibd.sc.samples$mild.v.other))
pval <- formatC(f7$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe/Moderate v Mild",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f7

# Disease severity: severe vs mild

table(ibd.sc.samples.2$SC.MHCI.Group,ibd.sc.samples.2$mild.v.severe)

f8 <- fisher.test(table(ibd.sc.samples.2$SC.MHCI.Group,ibd.sc.samples.2$mild.v.severe))
pval <- formatC(f8$p.value,digits=2,format="e")
update <- data.frame(Feature="Severe v Mild (subset)",Low.v.Other.Fisher=pval)
sc.fisher <- rbind(sc.fisher,update)

f8


##########################################################################
# CLINICAL DATA CORRELATIONS (2 CLUSTERS): IBD ONLY - WILCOX SCORE TESTS #
##########################################################################

# CLINICAL DATA CORRELATIONS (Average MHCI Score):
# TI/SC IBD samples only (ie everything except the controls)

# TI SAMPLES:
# TI Inflammation:

ibd.samples$TI.Inflammation.2 <- ifelse(ibd.samples$TI.numeric==1,"Y","N")

w1 <- wilcox_test(ibd.samples, Avg.MHCI.CpG~TI.numeric, p.adjust.method="BH")

w1

w1 <- w1 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-IBD_AvgMHCI_Inflammation_31Oct22.pdf")
ggplot(ibd.samples, aes(x=TI.Inflammation.2, y=Avg.MHCI.CpG, group=TI.Inflammation.2))+geom_boxplot(aes(fill=TI.Inflammation.2), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("TI Inflammation")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=TI.Inflammation.2), color="black",shape=21,size=2,alpha=0.9)+fillscale_inflammation2+stat_pvalue_manual(w1, label="p", tip.length=0.01)
dev.off()

ti.wilcox <- data.frame(Feature="Inflammation",AvgMHC.Wilcox=w1$p)

# Biologics:

w2 <- wilcox_test(ibd.samples, Avg.MHCI.CpG~Biologics.numeric, p.adjust.method="BH")

w2

w2 <- w2 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-IBD_AvgMHCI_Biologics_31Oct22.pdf")
ggplot(ibd.samples, aes(x=Biologics, y=Avg.MHCI.CpG, group=Biologics))+geom_boxplot(aes(fill=Biologics), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Biologics")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Biologics), color="black",shape=21,size=2,alpha=0.9)+fillscale_biologics+stat_pvalue_manual(w2, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Biologics",AvgMHC.Wilcox=w2$p)
ti.wilcox <- rbind(ti.wilcox,update)

# Surgery:

w3 <- wilcox_test(ibd.samples, Avg.MHCI.CpG~Surgery.numeric, p.adjust.method="BH")

w3

w3 <- w3 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-IBD_AvgMHCI_Surgery_31Oct22.pdf")
ggplot(ibd.samples, aes(x=Surgery, y=Avg.MHCI.CpG, group=Surgery))+geom_boxplot(aes(fill=Surgery), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Surgery")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Surgery), color="black",shape=21,size=2,alpha=0.9)+fillscale_surgery+stat_pvalue_manual(w3, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Surgery",AvgMHC.Wilcox=w3$p)
ti.wilcox <- rbind(ti.wilcox,update)


# AZA treatment:

w4 <- wilcox_test(ibd.samples, Avg.MHCI.CpG~AZA.numeric, p.adjust.method="BH")

w4

w4 <- w4 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-IBD_AvgMHCI_AZATreatment_31Oct22.pdf")
ggplot(ibd.samples, aes(x=AZA, y=Avg.MHCI.CpG, group=AZA))+geom_boxplot(aes(fill=AZA), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("AZA")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=AZA), color="black",shape=21,size=2,alpha=0.9)+fillscale_immunosuppressor+stat_pvalue_manual(w4, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="AZA Treatment",AvgMHC.Wilcox=w4$p)
ti.wilcox <- rbind(ti.wilcox,update)

# Perianal disease:

w5 <- wilcox_test(ibd.samples, Avg.MHCI.CpG~Perianal.numeric, p.adjust.method="BH")

w5

w5 <- w5 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-IBD_AvgMHCI_PerianalDisease_31Oct22.pdf")
ggplot(ibd.samples, aes(x=Perianal.Disease, y=Avg.MHCI.CpG, group=Perianal.Disease))+geom_boxplot(aes(fill=Perianal.Disease), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Perianal Disease")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Perianal.Disease), color="black",shape=21,size=2,alpha=0.9)+fillscale_perianal+stat_pvalue_manual(w5, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Perianal Disease",AvgMHC.Wilcox=w5$p)
ti.wilcox <- rbind(ti.wilcox,update)

# Disease severity: severe vs mild/moderate

w6 <- wilcox_test(ibd.samples, Avg.MHCI.CpG~severe.v.other, p.adjust.method="BH")

w6

w6 <- w6 %>% add_y_position()
ibd.samples$severe.v.other <- as.character(ibd.samples$severe.v.other)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-IBD_AvgMHCI_Severe-Other_31Oct22.pdf")
ggplot(ibd.samples, aes(x=severe.v.other, y=Avg.MHCI.CpG, group=severe.v.other))+geom_boxplot(aes(fill=severe.v.other), width=0.7)+theme_bw()+scale_x_discrete(limits=c("1","0"), labels=c("1"="Severe","0"="Mild/Moderate"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=severe.v.other), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity2+stat_pvalue_manual(w6, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe v Mild/Moderate",AvgMHC.Wilcox=w6$p)
ti.wilcox <- rbind(ti.wilcox,update)

# Disease severity: severe/moderate vs mild

w7 <- wilcox_test(ibd.samples, Avg.MHCI.CpG~mild.v.other, p.adjust.method="BH")

w7

w7 <- w7 %>% add_y_position()
ibd.samples$mild.v.other <- as.character(ibd.samples$mild.v.other)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-IBD_AvgMHCI_Mild-Other_31Oct22.pdf")
ggplot(ibd.samples, aes(x=mild.v.other, y=Avg.MHCI.CpG, group=mild.v.other))+geom_boxplot(aes(fill=mild.v.other), width=0.7)+theme_bw()+scale_x_discrete(limits=c("0","1"), labels=c("0"="Moderate/Severe","1"="Mild"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=mild.v.other), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity3+stat_pvalue_manual(w7, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe/Moderate v Mild",AvgMHC.Wilcox=w7$p)
ti.wilcox <- rbind(ti.wilcox,update)

# Disease severity: severe vs mild

w8 <- wilcox_test(ibd.samples.2, Avg.MHCI.CpG~mild.v.severe, p.adjust.method="BH")

w8

w8 <- w8 %>% add_y_position()
ibd.samples.2$mild.v.severe <- as.character(ibd.samples.2$mild.v.severe)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-IBD_AvgMHCI_Mild-Severe_31Oct22.pdf")
ggplot(ibd.samples.2, aes(x=mild.v.severe, y=Avg.MHCI.CpG, group=mild.v.severe))+geom_boxplot(aes(fill=mild.v.severe), width=0.7)+theme_bw()+scale_x_discrete(limits=c("1","0"), labels=c("1"="Severe","0"="Mild"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=mild.v.severe), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity2+stat_pvalue_manual(w8, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe v Mild (subset)",AvgMHC.Wilcox=w8$p)
ti.wilcox <- rbind(ti.wilcox,update)

ti.stats <- join(ti.wilcox, ti.fisher, type="left", match="all")
write.table(ti.stats, "MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-IBD_MHCI_Correlations_31Oct22.txt", row.names=F, sep="\t", quote=F)


# SC SAMPLES:
# SC Inflammation:

ibd.sc.samples$SC.Inflammation.2 <- ifelse(ibd.sc.samples$SC.numeric==1,"Y","N")

w1 <- wilcox_test(ibd.sc.samples, Avg.MHCI.CpG~SC.numeric, p.adjust.method="BH")

w1

w1 <- w1 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-IBD_AvgMHCI_Inflammation_31Oct22.pdf")
ggplot(ibd.sc.samples, aes(x=SC.Inflammation.2, y=Avg.MHCI.CpG, group=SC.Inflammation.2))+geom_boxplot(aes(fill=SC.Inflammation.2), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("SC Inflammation")+ylab("Average MHCI Methylation (SC)")+geom_jitter(aes(fill=SC.Inflammation.2), color="black",shape=21,size=2,alpha=0.9)+fillscale_inflammation2+stat_pvalue_manual(w1, label="p", tip.length=0.01)
dev.off()

sc.wilcox <- data.frame(Feature="Inflammation",AvgMHC.Wilcox=w1$p)

# Biologics:

w2 <- wilcox_test(ibd.sc.samples, Avg.MHCI.CpG~Biologics.numeric, p.adjust.method="BH")

w2

w2 <- w2 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-IBD_AvgMHCI_Biologics_31Oct22.pdf")
ggplot(ibd.sc.samples, aes(x=Biologics, y=Avg.MHCI.CpG, group=Biologics))+geom_boxplot(aes(fill=Biologics), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Biologics")+ylab("Average MHCI Methylation (SC)")+geom_jitter(aes(fill=Biologics), color="black",shape=21,size=2,alpha=0.9)+fillscale_biologics+stat_pvalue_manual(w2, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Biologics",AvgMHC.Wilcox=w2$p)
sc.wilcox <- rbind(sc.wilcox,update)

# Surgery:

w3 <- wilcox_test(ibd.sc.samples, Avg.MHCI.CpG~Surgery.numeric, p.adjust.method="BH")

w3

w3 <- w3 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-IBD_AvgMHCI_Surgery_31Oct22.pdf")
ggplot(ibd.sc.samples, aes(x=Surgery, y=Avg.MHCI.CpG, group=Surgery))+geom_boxplot(aes(fill=Surgery), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Surgery")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Surgery), color="black",shape=21,size=2,alpha=0.9)+fillscale_surgery+stat_pvalue_manual(w3, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Surgery",AvgMHC.Wilcox=w3$p)
sc.wilcox <- rbind(sc.wilcox,update)


# AZA treatment:

w4 <- wilcox_test(ibd.sc.samples, Avg.MHCI.CpG~AZA.numeric, p.adjust.method="BH")

w4

w4 <- w4 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-IBD_AvgMHCI_AZATreatment_31Oct22.pdf")
ggplot(ibd.sc.samples, aes(x=AZA, y=Avg.MHCI.CpG, group=AZA))+geom_boxplot(aes(fill=AZA), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("AZA")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=AZA), color="black",shape=21,size=2,alpha=0.9)+fillscale_immunosuppressor+stat_pvalue_manual(w4, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="AZA Treatment",AvgMHC.Wilcox=w4$p)
sc.wilcox <- rbind(sc.wilcox,update)

# Perianal disease:

w5 <- wilcox_test(ibd.sc.samples, Avg.MHCI.CpG~Perianal.numeric, p.adjust.method="BH")

w5

w5 <- w5 %>% add_y_position()

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-IBD_AvgMHCI_PerianalDisease_31Oct22.pdf")
ggplot(ibd.sc.samples, aes(x=Perianal.Disease, y=Avg.MHCI.CpG, group=Perianal.Disease))+geom_boxplot(aes(fill=Perianal.Disease), width=0.7)+theme_bw()+scale_x_discrete(limits=c("Y","N"), labels=c("N"="No","Y"="Yes"))+xlab("Perianal Disease")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Perianal.Disease), color="black",shape=21,size=2,alpha=0.9)+fillscale_perianal+stat_pvalue_manual(w5, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Perianal Disease",AvgMHC.Wilcox=w5$p)
sc.wilcox <- rbind(sc.wilcox,update)

# Disease severity: severe vs mild/moderate

w6 <- wilcox_test(ibd.sc.samples, Avg.MHCI.CpG~severe.v.other, p.adjust.method="BH")

w6

w6 <- w6 %>% add_y_position()
ibd.sc.samples$severe.v.other <- as.character(ibd.sc.samples$severe.v.other)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-IBD_AvgMHCI_Severe-Other_31Oct22.pdf")
ggplot(ibd.sc.samples, aes(x=severe.v.other, y=Avg.MHCI.CpG, group=severe.v.other))+geom_boxplot(aes(fill=severe.v.other), width=0.7)+theme_bw()+scale_x_discrete(limits=c("1","0"), labels=c("1"="Severe","0"="Mild/Moderate"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=severe.v.other), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity2+stat_pvalue_manual(w6, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe v Mild/Moderate",AvgMHC.Wilcox=w6$p)
sc.wilcox <- rbind(sc.wilcox,update)

# Disease severity: severe/moderate vs mild

w7 <- wilcox_test(ibd.sc.samples, Avg.MHCI.CpG~mild.v.other, p.adjust.method="BH")

w7

w7 <- w7 %>% add_y_position()
ibd.sc.samples$mild.v.other <- as.character(ibd.sc.samples$mild.v.other)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-IBD_AvgMHCI_Mild-Other_31Oct22.pdf")
ggplot(ibd.sc.samples, aes(x=mild.v.other, y=Avg.MHCI.CpG, group=mild.v.other))+geom_boxplot(aes(fill=mild.v.other), width=0.7)+theme_bw()+scale_x_discrete(limits=c("0","1"), labels=c("0"="Moderate/Severe","1"="Mild"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=mild.v.other), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity3+stat_pvalue_manual(w7, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe/Moderate v Mild",AvgMHC.Wilcox=w7$p)
sc.wilcox <- rbind(sc.wilcox,update)

# Disease severity: severe vs mild

w8 <- wilcox_test(ibd.sc.samples.2, Avg.MHCI.CpG~mild.v.severe, p.adjust.method="BH")

w8

w8 <- w8 %>% add_y_position()
ibd.sc.samples.2$mild.v.severe <- as.character(ibd.sc.samples.2$mild.v.severe)

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-IBD_AvgMHCI_Mild-Severe_31Oct22.pdf")
ggplot(ibd.sc.samples.2, aes(x=mild.v.severe, y=Avg.MHCI.CpG, group=mild.v.severe))+geom_boxplot(aes(fill=mild.v.severe), width=0.7)+theme_bw()+scale_x_discrete(limits=c("1","0"), labels=c("1"="Severe","0"="Mild"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=mild.v.severe), color="black",shape=21,size=2,alpha=0.9)+fillscale_severity2+stat_pvalue_manual(w8, label="p", tip.length=0.01)
dev.off()

update <- data.frame(Feature="Severe v Mild (subset)",AvgMHC.Wilcox=w8$p)
sc.wilcox <- rbind(sc.wilcox,update)

sc.stats <- join(sc.wilcox, sc.fisher, type="left", match="all")
write.table(sc.stats, "MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-IBD_MHCI_Correlations_31Oct22.txt", row.names=F, sep="\t", quote=F)


#################################################
# AVERAGE METHYLATION SCORE vs DISEASE SEVERITY #
#################################################

# Improved plots
# TI SAMPLES:

ti.avg.1 <- cd.samples[,c("array.id","severe.v.other","Avg.MHCI.CpG")]
ti.avg.1 <- melt(ti.avg.1, c("array.id","severe.v.other"), na.rm=T)
colnames(ti.avg.1) <- c("Sample","Severity","CpG","Avg.MHCI")

p1 <- ggplot(ti.avg.1, aes(x=Severity, y=Avg.MHCI, group=Severity))+geom_boxplot(aes(fill=Severity), width=0.7)
s1 <- ti.avg.1 %>% wilcox_test(Avg.MHCI~Severity) %>% adjust_pvalue(method="BH")
s1 <- s1 %>% add_y_position()

s1

s1$p.sci <- ifelse(s1$p.adj<=0.05, (formatC(s1$p.adj, digits=2, format="e")), "ns")

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-CD_AvgMHCI-v-DiseaseSeverity_31Oct22.pdf", height=6, width=5)
p1+theme_bw()+scale_x_discrete(limits=c("severe","other"), labels=c("severe"="severe","other"="mild/moderate"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (TI)")+geom_jitter(aes(fill=Severity),color="black",shape=21,size=2,alpha=0.9,width=0.2)+fillscale_2_severity+stat_pvalue_manual(s1, label="p.sci", tip.length=0.01)
dev.off()


# SC SAMPLES:

sc.avg.1 <- cd.sc.samples[,c("array.id","severe.v.other","Avg.MHCI.CpG")]
sc.avg.1 <- melt(sc.avg.1, c("array.id","severe.v.other"), na.rm=T)
colnames(sc.avg.1) <- c("Sample","Severity","CpG","Avg.MHCI")

p1 <- ggplot(sc.avg.1, aes(x=Severity, y=Avg.MHCI, group=Severity))+geom_boxplot(aes(fill=Severity), width=0.7)
s1 <- sc.avg.1 %>% wilcox_test(Avg.MHCI~Severity) %>% adjust_pvalue(method="BH")
s1 <- s1 %>% add_y_position()

s1

s1$p.sci <- ifelse(s1$p.adj<=0.05, (formatC(s1$p.adj, digits=2, format="e")), "ns")

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-CD_AvgMHCI-v-DiseaseSeverity_31Oct22.pdf", height=6, width=5)
p1+theme_bw()+scale_x_discrete(limits=c("severe","other"), labels=c("severe"="severe","other"="mild/moderate"))+xlab("Disease Severity")+ylab("Average MHCI Methylation (SC)")+geom_jitter(aes(fill=Severity),color="black",shape=21,size=2,alpha=0.9,width=0.2)+fillscale_2_severity+stat_pvalue_manual(s1, label="p.sci", tip.length=0.01)
dev.off()


#############################################
# TOP NLRC5 METHYLATION vs DISEASE SEVERITY #
#############################################

# TI SAMPLES:

ti.nlrc5 <- ti.samples[,c("array.id","diagnosis","severe.v.other","cg07839457","cg07862320")]
ti.nlrc5 <- melt(ti.nlrc5, c("array.id","diagnosis","severe.v.other"), na.rm=T)
colnames(ti.nlrc5) <- c("Sample","Diagnosis","Severity","CpG","Beta")
ti.nlrc5.2 <- subset(ti.nlrc5, Diagnosis!="Control")

p1 <- ggplot(ti.nlrc5.2, aes(x=Severity, y=Beta, group=Severity))+geom_boxplot(aes(fill=Severity), width=0.7, outlier.shape=NA)
s1 <- ti.nlrc5.2 %>% dplyr::group_by(CpG) %>% wilcox_test(Beta~Severity) %>% adjust_pvalue(method="BH")
s1 <- s1 %>% add_y_position()

s1

s1$p.sci <- ifelse(s1$p.adj<=0.05, (formatC(s1$p.adj, digits=2, format="e")), "ns")

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_TI-CD_NLRC5-v-DiseaseSeverity_31Oct22.pdf", height=6, width=11)
p1+facet_grid(cols=vars(CpG))+theme_bw()+scale_x_discrete(limits=c("severe","other"), labels=c("severe"="severe","other"="mild/moderate"))+xlab("Disease Severity")+ylab("Beta (TI)")+geom_jitter(aes(fill=Severity),color="black",shape=21,size=2,alpha=0.9,width=0.2)+fillscale_2_severity+stat_pvalue_manual(s1, label="p.sci", tip.length=0.01)
dev.off()


# SC SAMPLES:

sc.nlrc5 <- sc.samples[,c("array.id","diagnosis","severe.v.other","cg07839457","cg07862320")]
sc.nlrc5 <- melt(sc.nlrc5, c("array.id","diagnosis","severe.v.other"), na.rm=T)
colnames(sc.nlrc5) <- c("Sample","Diagnosis","Severity","CpG","Beta")
sc.nlrc5.2 <- subset(sc.nlrc5, Diagnosis!="Control")

p1 <- ggplot(sc.nlrc5.2, aes(x=Severity, y=Beta, group=Severity))+geom_boxplot(aes(fill=Severity), width=0.7, outlier.shape=NA)
s1 <- sc.nlrc5.2 %>% dplyr::group_by(CpG) %>% wilcox_test(Beta~Severity) %>% adjust_pvalue(method="BH")
s1 <- s1 %>% add_y_position()

s1

s1$p.sci <- ifelse(s1$p.adj<=0.05, (formatC(s1$p.adj, digits=2, format="e")), "ns")

pdf("MHCI_SCORE_CORRELATIONS/PaediatricOrganoid_SC-CD_NLRC5-v-DiseaseSeverity_31Oct22.pdf", height=6, width=11)
p1+facet_grid(cols=vars(CpG))+theme_bw()+scale_x_discrete(limits=c("severe","other"), labels=c("severe"="severe","other"="mild/moderate"))+xlab("Disease Severity")+ylab("Beta (SC)")+geom_jitter(aes(fill=Severity),color="black",shape=21,size=2,alpha=0.9,width=0.2)+fillscale_2_severity+stat_pvalue_manual(s1, label="p.sci", tip.length=0.01)
dev.off()


#############################
# EXPANDED TI CLUSTER PLOTS #
#############################

# Using pre-saved clustering input

# All TI samples:

ti.diagnosis <- ggplot()+geom_rect(aes(xmin=1:nrow(ti.samples), xmax=1:nrow(ti.samples)+1, ymin=0, ymax=1, fill=ti.samples$diagnosis[match(labels(dend.ti), ti.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_RE_diagnosis2+myTheme_meanmeth_label+guides(fill=guide_legend(title="Diagnosis: ", title.theme=element_text(size=22)))
ti.inflam <- ggplot()+geom_rect(aes(xmin=1:nrow(ti.samples), xmax=1:nrow(ti.samples)+1, ymin=0, ymax=1, fill=ti.samples$TI.Inflammation[match(labels(dend.ti), ti.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_inflammation+myTheme_meanmeth_label+guides(fill=guide_legend(title="TI Inflammation: ", title.theme=element_text(size=22)))
ti.biologics <- ggplot()+geom_rect(aes(xmin=1:nrow(ti.samples), xmax=1:nrow(ti.samples)+1, ymin=0, ymax=1, fill=ti.samples$Biologics[match(labels(dend.ti), ti.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_biologics+myTheme_meanmeth_label+guides(fill=guide_legend(title="Biologics: ", title.theme=element_text(size=22)))
ti.surgery <- ggplot()+geom_rect(aes(xmin=1:nrow(ti.samples), xmax=1:nrow(ti.samples)+1, ymin=0, ymax=1, fill=ti.samples$Surgery[match(labels(dend.ti), ti.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_surgery+myTheme_meanmeth_label+guides(fill=guide_legend(title="Surgery: ", title.theme=element_text(size=22)))
ti.aza <- ggplot()+geom_rect(aes(xmin=1:nrow(ti.samples), xmax=1:nrow(ti.samples)+1, ymin=0, ymax=1, fill=ti.samples$AZA[match(labels(dend.ti), ti.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_immunosuppressor+myTheme_meanmeth_label+guides(fill=guide_legend(title="AZA: ", title.theme=element_text(size=22)))
ti.perianal <- ggplot()+geom_rect(aes(xmin=1:nrow(ti.samples), xmax=1:nrow(ti.samples)+1, ymin=0, ymax=1, fill=ti.samples$Perianal.Disease[match(labels(dend.ti), ti.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_perianal+myTheme_meanmeth_label+guides(fill=guide_legend(title="Perianal Disease: ", title.theme=element_text(size=22)))
ti.severity <- ggplot()+geom_rect(aes(xmin=1:nrow(ti.samples), xmax=1:nrow(ti.samples)+1, ymin=0, ymax=1, fill=ti.samples$Disease.Severity[match(labels(dend.ti), ti.samples$array.id)]), color="black", alpha=0.8)+theme_bw()+fillscale_severity+myTheme_meanmeth_label+guides(fill=guide_legend(title="Disease Severity: ", title.theme=element_text(size=22)))
ti.meth <- ggplot(order.ti, aes(reorder(array.id,order), Avg.MHCI.CpG, group=1))+geom_point()+geom_line()+theme_bw()+myTheme_meanmeth_meth
ti.score <- ggplot()+geom_rect(aes(xmin=1:nrow(ti.samples), xmax=1:nrow(ti.samples)+1, ymin=0, ymax=1, fill=ti.samples$TI.MHCI.Group.2[match(labels(dend.ti), ti.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_3_cluster+myTheme_meanmeth_label+guides(fill=guide_legend(title="MHCI Group: ", title.theme=element_text(size=22)))

vp1 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.74, x=0.505)
vp2 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.67, x=0.505)
vp3 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.61, x=0.505)
vp4 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.55, x=0.505)
vp5 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.49, x=0.505)
vp6 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.43, x=0.505)
vp7 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.37, x=0.505)
vp8 <- viewport(height=unit(0.15,"npc"), width=unit(0.865,"npc"), just=c("center","top"), y=0.31, x=0.5)
vp9 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.15, x=0.505)

pdf("CLUSTER_PLOTS/PaediatricOrganoid_TI_ExpandedAverageMethylationClustering_31Oct22.pdf", width=27, height=20, onefile=F)
par(mfcol=c(1,1), mar=c(70,6,5,5)+0.1, oma=c(5,1,0,1))
plot.new()
myplclust(hc.ti, labels=ti.samples$case.no, lab.col=ti.samples$col.diagnosis, cex=1.2, main=NA, ylab=NA)
print(ti.diagnosis, vp=vp1)
grid.text("Diagnosis", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.722,"npc"), gp=gpar(fontsize=15))
print(ti.severity, vp=vp2)
grid.text("Disease\nSeverity", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.654,"npc"), gp=gpar(fontsize=15))
print(ti.inflam, vp=vp3)
grid.text("Inflammation", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.594,"npc"), gp=gpar(fontsize=15))
print(ti.biologics, vp=vp4)
grid.text("Biologics", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.532,"npc"), gp=gpar(fontsize=15))
print(ti.surgery, vp=vp5)
grid.text("Surgery", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.474,"npc"), gp=gpar(fontsize=15))
print(ti.aza, vp=vp6)
grid.text("AZA", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.414,"npc"), gp=gpar(fontsize=15))
print(ti.perianal, vp=vp7)
grid.text("Perianal\nDisease", just=c("right","center"), x=unit(0.064,"npc"), y=unit(0.355,"npc"), gp=gpar(fontsize=15))
print(ti.meth, vp=vp8)
grid.text("Average\nMHCI\nMethylation", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.25,"npc"), gp=gpar(fontsize=15))
print(ti.score, vp=vp9)
grid.text("MHCI Score", just=c("right","center"), x=unit(0.064,"npc"), y=unit(0.13,"npc"), gp=gpar(fontsize=15))
dev.off()

# CD TI samples only:

cd.samples <- subset(ti.samples, diagnosis=="CD")
cd.mhc <- ti.mhc[,cd.samples$array.id]

d.cd <- dist(t(cd.mhc))
hc.cd <- hclust(d.cd, method="complete")
dend.cd <- as.dendrogram(hc.cd)
order.cd <- cd.samples[match(labels(dend.cd),cd.samples$array.id),]
order.cd$order <- c(1:nrow(order.cd))
save(d.cd,hc.cd,dend.cd,order.cd,cd.samples,cd.mhc,file="CLUSTER_PLOTS/PaediatricOrganoid_TI-CD_ClusteringInput_31Oct22.RData")

ti.inflam <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.samples), xmax=1:nrow(cd.samples)+1, ymin=0, ymax=1, fill=cd.samples$TI.Inflammation[match(labels(dend.cd), cd.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_inflammation+myTheme_meanmeth_label+guides(fill=guide_legend(title="TI Inflammation: ", title.theme=element_text(size=22)))
ti.biologics <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.samples), xmax=1:nrow(cd.samples)+1, ymin=0, ymax=1, fill=cd.samples$Biologics[match(labels(dend.cd), cd.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_biologics+myTheme_meanmeth_label+guides(fill=guide_legend(title="Biologics: ", title.theme=element_text(size=22)))
ti.surgery <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.samples), xmax=1:nrow(cd.samples)+1, ymin=0, ymax=1, fill=cd.samples$Surgery[match(labels(dend.cd), cd.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_surgery+myTheme_meanmeth_label+guides(fill=guide_legend(title="Surgery: ", title.theme=element_text(size=22)))
ti.aza <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.samples), xmax=1:nrow(cd.samples)+1, ymin=0, ymax=1, fill=cd.samples$AZA[match(labels(dend.cd), cd.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_immunosuppressor+myTheme_meanmeth_label+guides(fill=guide_legend(title="AZA: ", title.theme=element_text(size=22)))
ti.perianal <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.samples), xmax=1:nrow(cd.samples)+1, ymin=0, ymax=1, fill=cd.samples$Perianal.Disease[match(labels(dend.cd), cd.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_perianal+myTheme_meanmeth_label+guides(fill=guide_legend(title="Perianal Disease: ", title.theme=element_text(size=22)))
ti.severity <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.samples), xmax=1:nrow(cd.samples)+1, ymin=0, ymax=1, fill=cd.samples$Disease.Severity[match(labels(dend.cd), cd.samples$array.id)]), color="black", alpha=0.8)+theme_bw()+fillscale_severity+myTheme_meanmeth_label+guides(fill=guide_legend(title="Disease Severity: ", title.theme=element_text(size=22)))
ti.meth <- ggplot(order.cd, aes(reorder(array.id,order), Avg.MHCI.CpG, group=1))+geom_point()+geom_line()+theme_bw()+myTheme_meanmeth_meth
ti.score <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.samples), xmax=1:nrow(cd.samples)+1, ymin=0, ymax=1, fill=cd.samples$TI.MHCI.Group.2[match(labels(dend.cd), cd.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_3_cluster+myTheme_meanmeth_label+guides(fill=guide_legend(title="MHCI Group: ", title.theme=element_text(size=22)))

vp2 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.67, x=0.505)
vp3 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.61, x=0.505)
vp4 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.55, x=0.505)
vp5 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.49, x=0.505)
vp6 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.43, x=0.505)
vp7 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.37, x=0.505)
vp9 <- viewport(height=unit(0.15,"npc"), width=unit(0.865,"npc"), just=c("center","top"), y=0.31, x=0.5)
vp10 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.15, x=0.505)

pdf("CLUSTER_PLOTS/PaediatricOrganoid_TI-CD_ExpandedAverageMethylationClustering_31Oct22.pdf", width=27, height=20, onefile=F)
par(mfcol=c(1,1), mar=c(70,6,5,5)+0.1, oma=c(5,1,0,1))
plot.new()
myplclust(hc.cd, labels=cd.samples$case.no, lab.col=cd.samples$col.diagnosis, cex=1.2, main=NA, ylab=NA)
print(ti.severity, vp=vp2)
grid.text("Disease\nSeverity", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.654,"npc"), gp=gpar(fontsize=15))
print(ti.inflam, vp=vp3)
grid.text("Inflammation", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.594,"npc"), gp=gpar(fontsize=15))
print(ti.biologics, vp=vp4)
grid.text("Biologics", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.532,"npc"), gp=gpar(fontsize=15))
print(ti.surgery, vp=vp5)
grid.text("Surgery", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.474,"npc"), gp=gpar(fontsize=15))
print(ti.aza, vp=vp6)
grid.text("AZA", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.414,"npc"), gp=gpar(fontsize=15))
print(ti.perianal, vp=vp7)
grid.text("Perianal\nDisease", just=c("right","center"), x=unit(0.064,"npc"), y=unit(0.355,"npc"), gp=gpar(fontsize=15))
print(ti.meth, vp=vp9)
grid.text("Average\nMHCI\nMethylation", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.25,"npc"), gp=gpar(fontsize=15))
print(ti.score, vp=vp10)
grid.text("MHCI Score", just=c("right","center"), x=unit(0.064,"npc"), y=unit(0.13,"npc"), gp=gpar(fontsize=15))
dev.off()


#############################
# EXPANDED SC CLUSTER PLOTS #
#############################

# Using pre-saved clustering input

# All SC samples:

sc.diagnosis <- ggplot()+geom_rect(aes(xmin=1:nrow(sc.samples), xmax=1:nrow(sc.samples)+1, ymin=0, ymax=1, fill=sc.samples$diagnosis[match(labels(dend.sc), sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_RE_diagnosis2+myTheme_meanmeth_label+guides(fill=guide_legend(title="Diagnosis: ", title.theme=element_text(size=22)))
sc.inflam <- ggplot()+geom_rect(aes(xmin=1:nrow(sc.samples), xmax=1:nrow(sc.samples)+1, ymin=0, ymax=1, fill=sc.samples$SC.Inflammation[match(labels(dend.sc), sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_inflammation+myTheme_meanmeth_label+guides(fill=guide_legend(title="SC Inflammation: ", title.theme=element_text(size=22)))
sc.biologics <- ggplot()+geom_rect(aes(xmin=1:nrow(sc.samples), xmax=1:nrow(sc.samples)+1, ymin=0, ymax=1, fill=sc.samples$Biologics[match(labels(dend.sc), sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_biologics+myTheme_meanmeth_label+guides(fill=guide_legend(title="Biologics: ", title.theme=element_text(size=22)))
sc.surgery <- ggplot()+geom_rect(aes(xmin=1:nrow(sc.samples), xmax=1:nrow(sc.samples)+1, ymin=0, ymax=1, fill=sc.samples$Surgery[match(labels(dend.sc), sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_surgery+myTheme_meanmeth_label+guides(fill=guide_legend(title="Surgery: ", title.theme=element_text(size=22)))
sc.aza <- ggplot()+geom_rect(aes(xmin=1:nrow(sc.samples), xmax=1:nrow(sc.samples)+1, ymin=0, ymax=1, fill=sc.samples$AZA[match(labels(dend.sc), sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_immunosuppressor+myTheme_meanmeth_label+guides(fill=guide_legend(title="AZA: ", title.theme=element_text(size=22)))
sc.perianal <- ggplot()+geom_rect(aes(xmin=1:nrow(sc.samples), xmax=1:nrow(sc.samples)+1, ymin=0, ymax=1, fill=sc.samples$Perianal.Disease[match(labels(dend.sc), sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_perianal+myTheme_meanmeth_label+guides(fill=guide_legend(title="Perianal Disease: ", title.theme=element_text(size=22)))
sc.severity <- ggplot()+geom_rect(aes(xmin=1:nrow(sc.samples), xmax=1:nrow(sc.samples)+1, ymin=0, ymax=1, fill=sc.samples$Disease.Severity[match(labels(dend.sc), sc.samples$array.id)]), color="black", alpha=0.8)+theme_bw()+fillscale_severity+myTheme_meanmeth_label+guides(fill=guide_legend(title="Disease Severity: ", title.theme=element_text(size=22)))
sc.meth <- ggplot(order.sc, aes(reorder(array.id,order), Avg.MHCI.CpG, group=1))+geom_point()+geom_line()+theme_bw()+myTheme_meanmeth_meth
sc.score <- ggplot()+geom_rect(aes(xmin=1:nrow(sc.samples), xmax=1:nrow(sc.samples)+1, ymin=0, ymax=1, fill=sc.samples$SC.MHCI.Group.2[match(labels(dend.sc), sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_3_cluster+myTheme_meanmeth_label+guides(fill=guide_legend(title="MHCI Group: ", title.theme=element_text(size=22)))

vp1 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.74, x=0.505)
vp2 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.67, x=0.505)
vp3 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.61, x=0.505)
vp4 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.55, x=0.505)
vp5 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.49, x=0.505)
vp6 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.43, x=0.505)
vp7 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.37, x=0.505)
vp8 <- viewport(height=unit(0.15,"npc"), width=unit(0.865,"npc"), just=c("center","top"), y=0.31, x=0.5)
vp9 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.15, x=0.505)

pdf("CLUSTER_PLOTS/PaediatricOrganoid_SC_ExpandedAverageMethylationClustering_31Oct22.pdf", width=27, height=20, onefile=F)
par(mfcol=c(1,1), mar=c(70,6,5,5)+0.1, oma=c(5,1,0,1))
plot.new()
myplclust(hc.sc, labels=sc.samples$case.no, lab.col=sc.samples$col.diagnosis, cex=1.2, main=NA, ylab=NA)
print(sc.diagnosis, vp=vp1)
grid.text("Diagnosis", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.722,"npc"), gp=gpar(fontsize=15))
print(sc.severity, vp=vp2)
grid.text("Disease\nSeverity", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.654,"npc"), gp=gpar(fontsize=15))
print(sc.inflam, vp=vp3)
grid.text("Inflammation", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.594,"npc"), gp=gpar(fontsize=15))
print(sc.biologics, vp=vp4)
grid.text("Biologics", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.532,"npc"), gp=gpar(fontsize=15))
print(sc.surgery, vp=vp5)
grid.text("Surgery", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.474,"npc"), gp=gpar(fontsize=15))
print(sc.aza, vp=vp6)
grid.text("AZA", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.414,"npc"), gp=gpar(fontsize=15))
print(sc.perianal, vp=vp7)
grid.text("Perianal\nDisease", just=c("right","center"), x=unit(0.064,"npc"), y=unit(0.355,"npc"), gp=gpar(fontsize=15))
print(sc.meth, vp=vp8)
grid.text("Average\nMHCI\nMethylation", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.25,"npc"), gp=gpar(fontsize=15))
print(sc.score, vp=vp9)
grid.text("MHCI Score", just=c("right","center"), x=unit(0.064,"npc"), y=unit(0.13,"npc"), gp=gpar(fontsize=15))
dev.off()


# CD SC samples only:

cd.sc.samples <- subset(sc.samples, diagnosis=="CD")
cd.sc.mhc <- sc.mhc[,cd.sc.samples$array.id]

d.cd.sc <- dist(t(cd.sc.mhc))
hc.cd.sc <- hclust(d.cd.sc, method="complete")
dend.cd.sc <- as.dendrogram(hc.cd.sc)
order.cd.sc <- cd.sc.samples[match(labels(dend.cd.sc),cd.sc.samples$array.id),]
order.cd.sc$order <- c(1:nrow(order.cd.sc))
save(d.cd.sc,hc.cd.sc,dend.cd.sc,order.cd.sc,cd.sc.samples,cd.sc.mhc,file="CLUSTER_PLOTS/PaediatricOrganoid_SC-CD_ClusteringInput_31Oct22.RData")

sc.inflam <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.sc.samples), xmax=1:nrow(cd.sc.samples)+1, ymin=0, ymax=1, fill=cd.sc.samples$SC.Inflammation[match(labels(dend.cd.sc), cd.sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_inflammation+myTheme_meanmeth_label+guides(fill=guide_legend(title="SC Inflammation: ", title.theme=element_text(size=22)))
sc.biologics <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.sc.samples), xmax=1:nrow(cd.sc.samples)+1, ymin=0, ymax=1, fill=cd.sc.samples$Biologics[match(labels(dend.cd.sc), cd.sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_biologics+myTheme_meanmeth_label+guides(fill=guide_legend(title="Biologics: ", title.theme=element_text(size=22)))
sc.surgery <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.sc.samples), xmax=1:nrow(cd.sc.samples)+1, ymin=0, ymax=1, fill=cd.sc.samples$Surgery[match(labels(dend.cd.sc), cd.sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_surgery+myTheme_meanmeth_label+guides(fill=guide_legend(title="Surgery: ", title.theme=element_text(size=22)))
sc.aza <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.sc.samples), xmax=1:nrow(cd.sc.samples)+1, ymin=0, ymax=1, fill=cd.sc.samples$AZA[match(labels(dend.cd.sc), cd.sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_immunosuppressor+myTheme_meanmeth_label+guides(fill=guide_legend(title="AZA: ", title.theme=element_text(size=22)))
sc.perianal <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.sc.samples), xmax=1:nrow(cd.sc.samples)+1, ymin=0, ymax=1, fill=cd.sc.samples$Perianal.Disease[match(labels(dend.cd.sc), cd.sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_perianal+myTheme_meanmeth_label+guides(fill=guide_legend(title="Perianal Disease: ", title.theme=element_text(size=22)))
sc.severity <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.sc.samples), xmax=1:nrow(cd.sc.samples)+1, ymin=0, ymax=1, fill=cd.sc.samples$Disease.Severity[match(labels(dend.cd.sc), cd.sc.samples$array.id)]), color="black", alpha=0.8)+theme_bw()+fillscale_severity+myTheme_meanmeth_label+guides(fill=guide_legend(title="Disease Severity: ", title.theme=element_text(size=22)))
sc.meth <- ggplot(order.cd.sc, aes(reorder(array.id,order), Avg.MHCI.CpG, group=1))+geom_point()+geom_line()+theme_bw()+myTheme_meanmeth_meth
sc.score <- ggplot()+geom_rect(aes(xmin=1:nrow(cd.sc.samples), xmax=1:nrow(cd.sc.samples)+1, ymin=0, ymax=1, fill=cd.sc.samples$SC.MHCI.Group.2[match(labels(dend.cd.sc), cd.sc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_3_cluster+myTheme_meanmeth_label+guides(fill=guide_legend(title="MHCI Group: ", title.theme=element_text(size=22)))

vp2 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.67, x=0.505)
vp3 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.61, x=0.505)
vp4 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.55, x=0.505)
vp5 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.49, x=0.505)
vp6 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.43, x=0.505)
vp7 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.37, x=0.505)
vp9 <- viewport(height=unit(0.15,"npc"), width=unit(0.865,"npc"), just=c("center","top"), y=0.31, x=0.5)
vp10 <- viewport(height=unit(0.065,"npc"), width=unit(0.93,"npc"), just=c("center","top"), y=0.15, x=0.505)

pdf("CLUSTER_PLOTS/PaediatricOrganoid_SC-CD_ExpandedAverageMethylationClustering_31Oct22.pdf", width=27, height=20, onefile=F)
par(mfcol=c(1,1), mar=c(70,6,5,5)+0.1, oma=c(5,1,0,1))
plot.new()
myplclust(hc.cd.sc, labels=cd.sc.samples$case.no, lab.col=cd.sc.samples$col.diagnosis, cex=1.2, main=NA, ylab=NA)
print(sc.severity, vp=vp2)
grid.text("Disease\nSeverity", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.654,"npc"), gp=gpar(fontsize=15))
print(sc.inflam, vp=vp3)
grid.text("Inflammation", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.594,"npc"), gp=gpar(fontsize=15))
print(sc.biologics, vp=vp4)
grid.text("Biologics", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.532,"npc"), gp=gpar(fontsize=15))
print(sc.surgery, vp=vp5)
grid.text("Surgery", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.474,"npc"), gp=gpar(fontsize=15))
print(sc.aza, vp=vp6)
grid.text("AZA", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.414,"npc"), gp=gpar(fontsize=15))
print(sc.perianal, vp=vp7)
grid.text("Perianal\nDisease", just=c("right","center"), x=unit(0.064,"npc"), y=unit(0.355,"npc"), gp=gpar(fontsize=15))
print(sc.meth, vp=vp9)
grid.text("Average\nMHCI\nMethylation", just=c("right","center"), x=unit(0.065,"npc"), y=unit(0.25,"npc"), gp=gpar(fontsize=15))
print(sc.score, vp=vp10)
grid.text("MHCI Score", just=c("right","center"), x=unit(0.064,"npc"), y=unit(0.13,"npc"), gp=gpar(fontsize=15))
dev.off()


##################################
# EXPANDED TI + SC CLUSTER PLOTS #
##################################

# SC/TI only colour scheme:

myColours <- c("darkolivegreen3","cornflowerblue")
colour_possibilities <- c("SC","TI")
names(myColours) <- colour_possibilities
fillscale_sampsite4 <- scale_fill_manual(name="Tissue", values = myColours, drop = T)
colscale_sampsite4 <- scale_colour_manual(name="Tissue", values = myColours, drop = T)


# CREATE A FULL TI+SC DATASET:

tisc.samples <- subset(samples, tissue!="DUO")
tisc.mhc <- beta.mhc[,tisc.samples$array.id]

# Add an inflammation column specific to the tissue the sample comes from:

tisc.samples$Inflammation <- NA
tisc.samples$Inflammation[tisc.samples$tissue=="SC"] <- tisc.samples$SC.Inflammation[tisc.samples$tissue=="SC"]
tisc.samples$Inflammation[tisc.samples$tissue=="TI"] <- tisc.samples$TI.Inflammation[tisc.samples$tissue=="TI"]

# Add MHCI group columns specific to the tissue the sample comes from:

tisc.samples$TI.MHCI.Group.2 <- tisc.samples$TI.MHCI.Group
tisc.samples$TI.MHCI.Group.2[tisc.samples$TI.Cluster==1] <- "Intermediate"
tisc.samples$TI.MHCI.Group.2[tisc.samples$TI.Cluster==2] <- "High"
tisc.samples$SC.MHCI.Group.2 <- tisc.samples$SC.MHCI.Group
tisc.samples$SC.MHCI.Group.2[tisc.samples$SC.Cluster==2] <- "Intermediate"
tisc.samples$SC.MHCI.Group.2[tisc.samples$SC.Cluster==1] <- "High"

tisc.samples$MHCI.Group.2 <- NA
tisc.samples$MHCI.Group.2[tisc.samples$tissue=="SC"] <- tisc.samples$SC.MHCI.Group.2[tisc.samples$tissue=="SC"]
tisc.samples$MHCI.Group.2[tisc.samples$tissue=="TI"] <- tisc.samples$TI.MHCI.Group.2[tisc.samples$tissue=="TI"]

# Make sure that there's a column with the diagnosis colour scheme:

tisc.samples$col.diagnosis <- as.factor(tisc.samples$diagnosis)
levels(tisc.samples$col.diagnosis) <- c("dodgerblue3","lightgrey","#636B4D","darkgoldenrod1")
tisc.samples$col.diagnosis <- as.character(tisc.samples$col.diagnosis)

# Recode the controls to "N" for biologics, surgery, AZA & perianal disease (currently mostly, but not all "NA"):

tisc.samples$Biologics[tisc.samples$diagnosis=="Control"] <- "N"
tisc.samples$Surgery[tisc.samples$diagnosis=="Control"] <- "N"
tisc.samples$AZA[tisc.samples$diagnosis=="Control"] <- "N"
tisc.samples$Perianal.Disease[tisc.samples$diagnosis=="Control"] <- "N"


# SUBSET TO JUST CD/CONTROL SAMPLES:

tisc.samples.2 <- subset(tisc.samples, diagnosis %in% c("Control","CD"))
tisc.mhc.2 <- tisc.mhc[,tisc.samples.2$array.id]


# EXPANDED CLUSTER PLOT FOR SC + TI CD/CONTROL SAMPLES ONLY:

d.2 <- dist(t(tisc.mhc.2))
hc.2 <- hclust(d.2, method="complete")
dend.2 <- as.dendrogram(hc.2)
order.2 <- tisc.samples.2[match(labels(dend.2),tisc.samples.2$array.id),]
order.2$order <- c(1:nrow(order.2))
save(tisc.samples.2,d.2,hc.2,dend.2,order.2,file="CLUSTER_PLOTS/PaediatricOrganoid_TISC_CD-Control_ClusteringInput_31Oct22.RData")

tissue <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples.2), xmax=1:nrow(tisc.samples.2)+1, ymin=0, ymax=1, fill=tisc.samples.2$tissue[match(labels(dend.2), tisc.samples.2$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_sampsite4+myTheme_meanmeth_label+guides(fill=guide_legend(title="Tissue: ", title.theme=element_text(size=22)))
diagnosis <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples.2), xmax=1:nrow(tisc.samples.2)+1, ymin=0, ymax=1, fill=tisc.samples.2$diagnosis[match(labels(dend.2), tisc.samples.2$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_RE_diagnosis2+myTheme_meanmeth_label+guides(fill=guide_legend(title="Diagnosis: ", title.theme=element_text(size=22)))
inflam <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples.2), xmax=1:nrow(tisc.samples.2)+1, ymin=0, ymax=1, fill=tisc.samples.2$Inflammation[match(labels(dend.2), tisc.samples.2$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_inflammation+myTheme_meanmeth_label+guides(fill=guide_legend(title="Inflammation: ", title.theme=element_text(size=22)))
biologics <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples.2), xmax=1:nrow(tisc.samples.2)+1, ymin=0, ymax=1, fill=tisc.samples.2$Biologics[match(labels(dend.2), tisc.samples.2$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_biologics+myTheme_meanmeth_label+guides(fill=guide_legend(title="Biologics: ", title.theme=element_text(size=22)))
surgery <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples.2), xmax=1:nrow(tisc.samples.2)+1, ymin=0, ymax=1, fill=tisc.samples.2$Surgery[match(labels(dend.2), tisc.samples.2$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_surgery+myTheme_meanmeth_label+guides(fill=guide_legend(title="Surgery: ", title.theme=element_text(size=22)))
aza <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples.2), xmax=1:nrow(tisc.samples.2)+1, ymin=0, ymax=1, fill=tisc.samples.2$AZA[match(labels(dend.2), tisc.samples.2$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_immunosuppressor+myTheme_meanmeth_label+guides(fill=guide_legend(title="AZA: ", title.theme=element_text(size=22)))
perianal <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples.2), xmax=1:nrow(tisc.samples.2)+1, ymin=0, ymax=1, fill=tisc.samples.2$Perianal.Disease[match(labels(dend.2), tisc.samples.2$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_perianal+myTheme_meanmeth_label+guides(fill=guide_legend(title="Perianal Disease: ", title.theme=element_text(size=22)))
severity <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples.2), xmax=1:nrow(tisc.samples.2)+1, ymin=0, ymax=1, fill=tisc.samples.2$Disease.Severity[match(labels(dend.2), tisc.samples.2$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_severity+myTheme_meanmeth_label+guides(fill=guide_legend(title="Disease Severity: ", title.theme=element_text(size=22)))
meth <- ggplot(order.2, aes(reorder(array.id,order), Avg.MHCI.CpG, group=1))+geom_point()+geom_line()+theme_bw()+myTheme_meanmeth_meth
score <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples.2), xmax=1:nrow(tisc.samples.2)+1, ymin=0, ymax=1, fill=tisc.samples.2$MHCI.Group.2[match(labels(dend.2), tisc.samples.2$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_3_cluster+myTheme_meanmeth_label+guides(fill=guide_legend(title="MHCI Group: ", title.theme=element_text(size=22)))

vp1 <- viewport(height=unit(0.065,"npc"), width=unit(0.955,"npc"), just=c("center","top"), y=0.74, x=0.502)
vp2 <- viewport(height=unit(0.065,"npc"), width=unit(0.955,"npc"), just=c("center","top"), y=0.67, x=0.502)
vp3 <- viewport(height=unit(0.065,"npc"), width=unit(0.955,"npc"), just=c("center","top"), y=0.61, x=0.502)
vp4 <- viewport(height=unit(0.065,"npc"), width=unit(0.955,"npc"), just=c("center","top"), y=0.55, x=0.502)
vp5 <- viewport(height=unit(0.065,"npc"), width=unit(0.955,"npc"), just=c("center","top"), y=0.49, x=0.502)
vp6 <- viewport(height=unit(0.065,"npc"), width=unit(0.955,"npc"), just=c("center","top"), y=0.43, x=0.502)
vp7 <- viewport(height=unit(0.065,"npc"), width=unit(0.955,"npc"), just=c("center","top"), y=0.37, x=0.502)
vp8 <- viewport(height=unit(0.065,"npc"), width=unit(0.955,"npc"), just=c("center","top"), y=0.31, x=0.502)
vp9 <- viewport(height=unit(0.15,"npc"), width=unit(0.878,"npc"), just=c("center","top"), y=0.25, x=0.498)
vp10 <- viewport(height=unit(0.065,"npc"), width=unit(0.955,"npc"), just=c("center","top"), y=0.09, x=0.502)

pdf("CLUSTER_PLOTS/PaediatricOrganoid_TISC_CD-Control_ExpandedAverageMethylationClustering_31Oct22.pdf", width=35, height=20, onefile=F)
par(mfcol=c(1,1), mar=c(70,6,5,5)+0.1, oma=c(5,1,0,1))
plot.new()
myplclust(hc.2, labels=tisc.samples.2$case.no, lab.col=tisc.samples.2$col.diagnosis, cex=1.2, main=NA, ylab=NA)
print(tissue, vp=vp1)
grid.text("Tissue", just=c("right","center"), x=unit(0.06,"npc"), y=unit(0.722,"npc"), gp=gpar(fontsize=15))
print(diagnosis, vp=vp2)
grid.text("Diagnosis", just=c("right","center"), x=unit(0.06,"npc"), y=unit(0.654,"npc"), gp=gpar(fontsize=15))
print(severity, vp=vp3)
grid.text("Disease\nSeverity", just=c("right","center"), x=unit(0.06,"npc"), y=unit(0.594,"npc"), gp=gpar(fontsize=15))
print(inflam, vp=vp4)
grid.text("Inflammation", just=c("right","center"), x=unit(0.06,"npc"), y=unit(0.532,"npc"), gp=gpar(fontsize=15))
print(biologics, vp=vp5)
grid.text("Biologics", just=c("right","center"), x=unit(0.06,"npc"), y=unit(0.474,"npc"), gp=gpar(fontsize=15))
print(surgery, vp=vp6)
grid.text("Surgery", just=c("right","center"), x=unit(0.06,"npc"), y=unit(0.414,"npc"), gp=gpar(fontsize=15))
print(aza, vp=vp7)
grid.text("AZA", just=c("right","center"), x=unit(0.06,"npc"), y=unit(0.355,"npc"), gp=gpar(fontsize=15))
print(perianal, vp=vp8)
grid.text("Perianal\nDisease", just=c("right","center"), x=unit(0.06,"npc"), y=unit(0.295,"npc"), gp=gpar(fontsize=15))
print(meth, vp=vp9)
grid.text("Average\nMHCI\nMethylation", just=c("right","center"), x=unit(0.06,"npc"), y=unit(0.19,"npc"), gp=gpar(fontsize=15))
print(score, vp=vp10)
grid.text("MHCI Score", just=c("right","center"), x=unit(0.06,"npc"), y=unit(0.07,"npc"), gp=gpar(fontsize=15))
dev.off()


# EXPANDED CLUSTER PLOT FOR ALL SC + TI SAMPLES:

d.all <- dist(t(tisc.mhc))
hc.all <- hclust(d.all, method="complete")
dend.all <- as.dendrogram(hc.all)
order.all <- tisc.samples[match(labels(dend.all),tisc.samples$array.id),]
order.all$order <- c(1:nrow(order.all))
save(tisc.samples,d.all,hc.all,dend.all,order.all,file="CLUSTER_PLOTS/PaediatricOrganoid_TISC_ClusteringInput_31Oct22.RData")

tissue <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples), xmax=1:nrow(tisc.samples)+1, ymin=0, ymax=1, fill=tisc.samples$tissue[match(labels(dend.all), tisc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_sampsite4+myTheme_meanmeth_label+guides(fill=guide_legend(title="Tissue: ", title.theme=element_text(size=22)))
diagnosis <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples), xmax=1:nrow(tisc.samples)+1, ymin=0, ymax=1, fill=tisc.samples$diagnosis[match(labels(dend.all), tisc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_RE_diagnosis2+myTheme_meanmeth_label+guides(fill=guide_legend(title="Diagnosis: ", title.theme=element_text(size=22)))
inflam <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples), xmax=1:nrow(tisc.samples)+1, ymin=0, ymax=1, fill=tisc.samples$Inflammation[match(labels(dend.all), tisc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_inflammation+myTheme_meanmeth_label+guides(fill=guide_legend(title="Inflammation: ", title.theme=element_text(size=22)))
biologics <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples), xmax=1:nrow(tisc.samples)+1, ymin=0, ymax=1, fill=tisc.samples$Biologics[match(labels(dend.all), tisc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_biologics+myTheme_meanmeth_label+guides(fill=guide_legend(title="Biologics: ", title.theme=element_text(size=22)))
surgery <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples), xmax=1:nrow(tisc.samples)+1, ymin=0, ymax=1, fill=tisc.samples$Surgery[match(labels(dend.all), tisc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_surgery+myTheme_meanmeth_label+guides(fill=guide_legend(title="Surgery: ", title.theme=element_text(size=22)))
aza <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples), xmax=1:nrow(tisc.samples)+1, ymin=0, ymax=1, fill=tisc.samples$AZA[match(labels(dend.all), tisc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_immunosuppressor+myTheme_meanmeth_label+guides(fill=guide_legend(title="AZA: ", title.theme=element_text(size=22)))
perianal <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples), xmax=1:nrow(tisc.samples)+1, ymin=0, ymax=1, fill=tisc.samples$Perianal.Disease[match(labels(dend.all), tisc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_perianal+myTheme_meanmeth_label+guides(fill=guide_legend(title="Perianal Disease: ", title.theme=element_text(size=22)))
severity <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples), xmax=1:nrow(tisc.samples)+1, ymin=0, ymax=1, fill=tisc.samples$Disease.Severity[match(labels(dend.all), tisc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_severity+myTheme_meanmeth_label+guides(fill=guide_legend(title="Disease Severity: ", title.theme=element_text(size=22)))
meth <- ggplot(order.all, aes(reorder(array.id,order), Avg.MHCI.CpG, group=1))+geom_point()+geom_line()+theme_bw()+myTheme_meanmeth_meth
score <- ggplot()+geom_rect(aes(xmin=1:nrow(tisc.samples), xmax=1:nrow(tisc.samples)+1, ymin=0, ymax=1, fill=tisc.samples$MHCI.Group.2[match(labels(dend.all), tisc.samples$array.id)]), color="black", alpha=0.6)+theme_bw()+fillscale_3_cluster+myTheme_meanmeth_label+guides(fill=guide_legend(title="MHCI Group: ", title.theme=element_text(size=22)))

vp1 <- viewport(height=unit(0.065,"npc"), width=unit(0.973,"npc"), just=c("center","top"), y=0.74, x=0.502)
vp2 <- viewport(height=unit(0.065,"npc"), width=unit(0.973,"npc"), just=c("center","top"), y=0.67, x=0.502)
vp3 <- viewport(height=unit(0.065,"npc"), width=unit(0.973,"npc"), just=c("center","top"), y=0.61, x=0.502)
vp4 <- viewport(height=unit(0.065,"npc"), width=unit(0.973,"npc"), just=c("center","top"), y=0.55, x=0.502)
vp5 <- viewport(height=unit(0.065,"npc"), width=unit(0.973,"npc"), just=c("center","top"), y=0.49, x=0.502)
vp6 <- viewport(height=unit(0.065,"npc"), width=unit(0.973,"npc"), just=c("center","top"), y=0.43, x=0.502)
vp7 <- viewport(height=unit(0.065,"npc"), width=unit(0.973,"npc"), just=c("center","top"), y=0.37, x=0.502)
vp8 <- viewport(height=unit(0.065,"npc"), width=unit(0.973,"npc"), just=c("center","top"), y=0.31, x=0.502)
vp9 <- viewport(height=unit(0.15,"npc"), width=unit(0.89,"npc"), just=c("center","top"), y=0.25, x=0.498)
vp10 <- viewport(height=unit(0.065,"npc"), width=unit(0.973,"npc"), just=c("center","top"), y=0.09, x=0.502)

pdf("CLUSTER_PLOTS/PaediatricOrganoid_TISC_ExpandedAverageMethylationClustering_31Oct22.pdf", width=50, height=20, onefile=F)
par(mfcol=c(1,1), mar=c(70,6,5,5)+0.1, oma=c(5,1,0,1))
plot.new()
myplclust(hc.all, labels=tisc.samples$case.no, lab.col=tisc.samples$col.diagnosis, cex=1.2, main=NA, ylab=NA)
print(tissue, vp=vp1)
grid.text("Tissue", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.722,"npc"), gp=gpar(fontsize=15))
print(diagnosis, vp=vp2)
grid.text("Diagnosis", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.654,"npc"), gp=gpar(fontsize=15))
print(severity, vp=vp3)
grid.text("Disease\nSeverity", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.594,"npc"), gp=gpar(fontsize=15))
print(inflam, vp=vp4)
grid.text("Inflammation", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.532,"npc"), gp=gpar(fontsize=15))
print(biologics, vp=vp5)
grid.text("Biologics", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.474,"npc"), gp=gpar(fontsize=15))
print(surgery, vp=vp6)
grid.text("Surgery", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.414,"npc"), gp=gpar(fontsize=15))
print(aza, vp=vp7)
grid.text("AZA", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.355,"npc"), gp=gpar(fontsize=15))
print(perianal, vp=vp8)
grid.text("Perianal\nDisease", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.295,"npc"), gp=gpar(fontsize=15))
print(meth, vp=vp9)
grid.text("Average\nMHCI\nMethylation", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.19,"npc"), gp=gpar(fontsize=15))
print(score, vp=vp10)
grid.text("MHCI Score", just=c("right","center"), x=unit(0.05,"npc"), y=unit(0.07,"npc"), gp=gpar(fontsize=15))
dev.off()


################
# SESSION INFO #
################

sessionInfo()
rm(list=ls())

