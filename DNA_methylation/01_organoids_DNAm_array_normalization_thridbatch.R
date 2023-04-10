#'---
#'title: DNAm preprocessing and normalization
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages(library(minfi))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(RColorBrewer))
suppressMessages(library(here))

suppressMessages(library(dplyr))
suppressMessages(library(lmtest))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools))
suppressMessages(library(rafalib))
suppressMessages(library(cowplot))



options(stringsAsFactors = FALSE)

subsample<-50

#' ### Load Functions
source(here("general_functions","00_pretty_plots.R"))
suppressMessages(source(here("general_functions","00_Heat_scree_plot_generic.R")))


#' ### Load Data
path<-"DNAm/data/batch_three_raw"
#' epic.organoid<-read.csv(here(path, "METH_K_Nayak_20038_SAMPLESHEET.csv"), skip=7)
#' epic.organoid$array.id<-paste(epic.organoid$Sentrix_ID, epic.organoid$Sentrix_Position, sep="_")
#' 
#' 
#' #'### add gender
#' load(here("DNAm/data", "AllEPICOrganoid_FullMetadata_23Sep22.RData"))
#' sampleInfo<-epic.samples[which(epic.samples$array.id%in%epic.organoid$array.id),]
#' epic.organoid<-merge(epic.organoid,sampleInfo, by="array.id")
#' 
#' #' ### Normalize DNAm Arrays
#' here(path)
#' epic.organoid$array.id.path <- file.path(here(path,"RAW_iDAT_FILES_20038_17JUN22",epic.organoid$array.id))
#' 
#' epic.organoid$individual<-sapply(1:nrow(epic.organoid), function(x) strsplit(epic.organoid$Sample_Name[x],"_")[[1]][1])
#' epic.organoid$condition<-as.factor(epic.organoid$treatment)
#' epic.organoid$condition<-as.character(epic.organoid$condition)
#' epic.organoid$comparison<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$condition[x]%in%c("IFNg","LPS","AZA")){"UT"}else{"NT"})
#' epic.organoid$Segment<-epic.organoid$tissue
#' epic.organoid$Gender<-epic.organoid$sex
#' 
#' 
#' # multiple DMAP files common with epic so need to force https://support.bioconductor.org/p/97773/
#' rgset_organoid <- read.metharray(epic.organoid$array.id.path, verbose = FALSE,force=TRUE)
#' 
#' # Background and dye bias correction with noob thhrough funnorm implemented in minfi
#' #http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
#' MSet.illumina <- preprocessFunnorm(rgset_organoid, sex=epic.organoid$Gender)
#' organoid_beta<-getBeta(MSet.illumina)
#' 
#' print(paste("Samples available: ",ncol(organoid_beta),"; Probes available: ",nrow(organoid_beta),sep=""))
#' 
#' save(rgset_organoid, MSet.illumina,epic.organoid, file=here("DNAm/data/batch_three_raw","thridbatch_raw.RData"))
#' 
#' 



load(here("DNAm/data/batch_three_raw","thridbatch_raw.RData"))
organoid_beta<-getBeta(MSet.illumina)

#' ### Detection p values across all probes for each sample
avg_detPval <- colMeans(detectionP(rgset_organoid))
epic.organoid$det_pval<-avg_detPval

ggplot(epic.organoid)+geom_boxplot(aes(as.factor(sentrix.id), det_pval, fill=as.factor(sentrix.id)), outlier.shape = NA)+
  geom_point(aes(as.factor(sentrix.id), det_pval, group=Sample_Name, fill=as.factor(sentrix.id)), shape=21, color="black",
             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE)

ggsave(here("DNAm/figs","thirdbatch_detection_pvalue_organoids.pdf"), width=6, height=5)
ggsave(here("DNAm/figs/jpeg","thirdbatch_detection_pvalue_organoids.jpeg"), width=6, height=5)




#' Beta distribution before and after normalization

# extract raw beta values for plotting
beta_raw<-getBeta(rgset_organoid)
identical(colnames(beta_raw),epic.organoid$array.id)

set.seed(1)
sample_sub<-sample(1:nrow(organoid_beta), nrow(organoid_beta)/subsample)
Beta_melted<- melt(organoid_beta[sample_sub,])
Beta_melted_raw<- melt(beta_raw[sample_sub,])

# Remove NAs before plotting (otherwise get many non-inifnite warnings)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
Beta_Plot_raw<-Beta_melted_raw[which(!(is.na(Beta_melted_raw$value))),]

# Add meta data
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")
colnames(Beta_Plot_raw)<-c("CpG","ID","Beta")
Beta_Plot_raw<-merge(Beta_Plot_raw,epic.organoid, by.x="ID", by.y="array.id")

beta_dis_EPIC_raw<-ggplot(Beta_Plot_raw, aes(Beta, group=as.character(ID), color=as.character(Segment)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")

beta_dis_EPIC<-ggplot(Beta_Plot, aes(Beta, group=as.character(ID), color=as.character(Segment)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")


grid.arrange(beta_dis_EPIC_raw,beta_dis_EPIC)

ggsave(here("DNAm/figs","thirdbatch_beta_distribution_organoid.pdf"),grid.arrange(beta_dis_EPIC_raw,beta_dis_EPIC),  w=10, h=5)
ggsave(here("DNAm/figs/jpeg","thirdbatch_beta_distribution_organoid.jpeg"),grid.arrange(beta_dis_EPIC_raw,beta_dis_EPIC), w=10, h=5)


ggplot(Beta_Plot, aes(Beta, group=as.character(ID), color=as.character(AZA)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")

ggsave(here("DNAm/figs","thirdbatch_beta_distribution_organoid_AZA.pdf"), w=5, h=5)
ggsave(here("DNAm/figs/jpeg","thirdbatch_beta_distribution_organoid_AZA.jpeg"),grid.arrange(beta_dis_EPIC_raw,beta_dis_EPIC), w=5, h=5)





#'#### Confirm individuals ID with SNPs probes and clustering by DNAm

# remove rows with NAs
Betas_cluster<-organoid_beta[complete.cases(organoid_beta),]

d <- dist(t(Betas_cluster))
hc <- hclust(d, method = "complete") #single, complete, average, ward
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(epic.organoid$Segment), cex=1.5)

pdf(here("DNAm/figs","thirdbatch_cluster_wholeEPIC_organoid.pdf"), width=30)
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(epic.organoid$Segment), cex=1.5)
dev.off()


#' #### Genotyping Probes
SNPs <- getSnpBeta(rgset_organoid)
SNPs<-SNPs[complete.cases(SNPs),]# 65 cause one was all NA

SNPs<-SNPs[,which(colnames(SNPs)%in%epic.organoid$array.id)]
identical(colnames(SNPs),epic.organoid$array.id)

d <- dist(t(SNPs))
hc <- hclust(d, method = "complete") #single, complete, average, ward
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(as.character(epic.organoid$individual)), cex=1.5)

pdf(here("DNAm/figs","thirdbatch_cluster_snps_EPIC_organoid.pdf"), width=30)
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(as.character(epic.organoid$individual)), cex=1.5)
dev.off()



#' #### Sex clustering
#' Using the cg ID to chromosome annotation from illumina
#' https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html
anno_EPIC<-read.csv(here("../data", "MethylationEPIC_v-1-0_B4.csv"), skip=7)

organoid_beta<-organoid_beta[which(rownames(organoid_beta)%in%anno_EPIC$IlmnID),]
anno_EPIC<-anno_EPIC[match(rownames(organoid_beta),anno_EPIC$IlmnID),]
identical(rownames(organoid_beta),anno_EPIC$IlmnID)

organoid_beta_sex<-organoid_beta[which(anno_EPIC$CHR%in%c('X','Y')),]

d <- dist(t(organoid_beta_sex))
hc <- hclust(d, method = "complete") #single, complete, average, ward
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(epic.organoid$Gender), cex=1.5)

pdf(here("DNAm/figs","thirdbatch_cluster_sex_EPIC_organoid.pdf"), width=30)
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(epic.organoid$Gender), cex=1.5)
dev.off()

#'##' Another cehck for sex mixups
epic.organoid$sex_mix<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$individual[x] %in% c("T036","T279","T203","T193","T238","T202","T196","T192","T362")){
  epic.organoid$Sample_Name[x]}else{epic.organoid$Gender[x]})

organoid_beta_sex <- organoid_beta[anno_EPIC$CHR%in%c("X", "Y"), ]
save(organoid_beta_sex, epic.organoid, file=here("DNAm/data/batch_three_raw","thirdbatch_betas_normalized_sexonly.RData"))

load(file=here("DNAm/data/batch_three_raw","thirdbatch_betas_normalized_sexonly.RData"))


# beta plot y chr
y_cpg<-anno_EPIC$IlmnID[which(anno_EPIC$CHR%in%c('Y'))]
ibd_Y<-organoid_beta_sex[which(rownames(organoid_beta_sex)%in%y_cpg),]#537
Beta_melted<- melt(ibd_Y)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")

all_sample<-ggplot(Beta_Plot, aes(Beta, group=Sample_Name, color=sex_mix))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")+
  scale_color_manual(values=c("#9ecae1","#fc9272", brewer.pal(length(which(!(epic.organoid$sex_mix%in%c("F","M")))), "Set3") ), name="Sex")
all_sample

Beta_Plot_mean<-Beta_Plot[which((Beta_Plot$sex_mix%in%c("M","F"))),]
Beta_Plot_mean$sex_mix<-NULL
sample_split<-ggplot()+
  geom_density(aes(Beta, group=sex_mix),Beta_Plot[which(!(Beta_Plot$sex_mix%in%c("M","F"))),])+
  geom_density(aes(Beta, group=Gender, color=Gender),Beta_Plot_mean)+
  theme_bw()+xlab("DNAm Beta Value")+facet_wrap(~sex_mix)+scale_color_manual(values=c("#9ecae1","#fc9272"))
sample_split

ggsave(here("DNAm/figs","thirdbatch_beta_distribution_EPIC_Y.pdf"),grid.arrange(all_sample,sample_split), w=10, h=15)
ggsave(here("DNAm/figs/jpeg","thirdbatch_beta_distribution_EPIC_Y.jpeg"),grid.arrange(all_sample,sample_split), w=10, h=15)

# beta plot x chr
x_cpg<-anno_EPIC$IlmnID[which(anno_EPIC$CHR%in%c('X'))]
ibd_X<-organoid_beta_sex[which(rownames(organoid_beta_sex)%in%x_cpg),]#19090
Beta_melted<- melt(ibd_X)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")

all_sample<-ggplot(Beta_Plot, aes(Beta, group=Sample_Name, color=sex_mix))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")+
  scale_color_manual(values=c("#9ecae1","#fc9272", brewer.pal(length(which(!(epic.organoid$sex_mix%in%c("F","M")))), "Set3") ), name="Sex")
all_sample

Beta_Plot_mean<-Beta_Plot[which((Beta_Plot$sex_mix%in%c("M","F"))),]
Beta_Plot_mean$sex_mix<-NULL
sample_split<-ggplot()+
  geom_density(aes(Beta, group=sex_mix),Beta_Plot[which(!(Beta_Plot$sex_mix%in%c("M","F"))),])+
  geom_density(aes(Beta, group=Gender, color=Gender),Beta_Plot_mean)+
  theme_bw()+xlab("DNAm Beta Value")+facet_wrap(~sex_mix)+scale_color_manual(values=c("#9ecae1","#fc9272"))
sample_split

ggsave(here("DNAm/figs","thirdbatch_beta_distribution_EPIC_X.pdf"),grid.arrange(all_sample,sample_split), w=10, h=15)
ggsave(here("DNAm/figs/jpeg","thirdbatch_beta_distribution_EPIC_X.jpeg"),grid.arrange(all_sample,sample_split), w=10, h=15)




#' #### Remove samples which do not cluster correctly
#' none to remove here though

#' ### Probe Filtering
# SNP probes should already be removed
organoid_beta <- organoid_beta[!grepl("rs",rownames(organoid_beta)), ]
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))

#' #### Sex Chromosomes
anno_EPIC<-anno_EPIC[anno_EPIC$IlmnID%in%rownames(organoid_beta),]
identical(rownames(organoid_beta),anno_EPIC$IlmnID)
organoid_beta <- organoid_beta[!anno_EPIC$CHR%in%c("X", "Y"), ]

filt_sex<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))


#' #### Cross-hybridizing probes and polymorphic probes.
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1
#' "43,254 cross-reactive probes with ≥ 47 bp homology with an off-target site, of which 15,782 (36.5 %) are new to the EPIC platform"
#' They include this annotated list in their supplement.
#' wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM2_ESM.csv
cross_reactive<-read.csv(here("../data","13059_2016_1066_MOESM2_ESM.csv"), stringsAsFactors = F)
organoid_beta<-organoid_beta[which(!(rownames(organoid_beta)%in%cross_reactive$PROBE)),]

filt_cross<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))


#'For polymorphic probes I will The Pidsley annotation aswell for "Probes overlapping genetic variants at targeted CpG sites." and "Probes overlapping genetic variants at single base extension sites for Infinium Type I probes" but NOT "Probes with genetic variants overlapping the body of the probe: 48 base pairs for Infinium Type I probes and 49 base pairs for Infinium Type II probes."

#wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM4_ESM.csv
#wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM5_ESM.csv

polymorphic<-read.csv(here("../data","13059_2016_1066_MOESM4_ESM.csv"),stringsAsFactors = F)
print(paste("Filtering ",length(unique(polymorphic$PROBE))," polymorphic probes (genetic variants at targeted CpG sites).", sep=""))

baseext<-read.csv(here("../data","13059_2016_1066_MOESM5_ESM.csv"),stringsAsFactors = F)
print(paste("Filtering ",length(unique(baseext$PROBE))," polymorphic probes (single base extension sites for Infinium Type I probes).", sep=""))

organoid_beta<-organoid_beta[which(!(rownames(organoid_beta)%in%c(polymorphic$PROBE, baseext$PROBE))),]

filt_poly<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))



#' #### Probe filtering based on detection pvalue and detection over background (NA)

#' Remove probes with high NA count
na_count_probe <-sapply(1:nrow(organoid_beta), function(y) length(which(is.na(organoid_beta[y,]))))
na_count_probe_good<-which(na_count_probe<(ncol(organoid_beta)*0.05))
organoid_beta<-organoid_beta[na_count_probe_good,]

filt_bead<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))


#' Remove probes with high detection p value across samples, and any samples with many high detection p value probes
detP <- detectionP(rgset_organoid)
detP<-detP[,which(colnames(detP)%in%epic.organoid$array.id)]
identical(colnames(detP),epic.organoid$array.id)

failed <- detP>0.05
bad_det_p<-names(which(rowMeans(failed)>0.01))
bad_det_psamp<-names(which(colMeans(failed)>0.01))

organoid_beta<-organoid_beta[which(!(rownames(organoid_beta)%in%bad_det_p)),]
organoid_beta<-organoid_beta[,which(!(colnames(organoid_beta)%in%bad_det_psamp))]
identical(colnames(organoid_beta), as.character(epic.organoid$array.id))

filt_detp<-nrow(organoid_beta)
print(paste("Samples available: ",ncol(organoid_beta),"\nProbes available: ",nrow(organoid_beta),sep=""))



#' #### Probe attrition plot
df<-data.frame(sample_num_remaining=c(866238,865918,865859,filt_sex,filt_cross,filt_poly,filt_bead,filt_detp),
               filter=c("EPIC Probe Number","Missing Annotation Data","Removal of SNP Probes",
                        "Removal of X and Y chromosome probes","Removal of Cross Reactive Probes",
                        "Removal of Polymorphic Probes", "Removal of Probes with Beadcount <3\nin 5 % of Samples",
                        "Removal of Probes with 1 % of samples\nwith a detection p-value greater than 0.05"))
df$sample_num_lost<-c(0,sapply(2:nrow(df), function(x) df$sample_num_remaining[x-1]-df$sample_num_remaining[x]))

df$filter<-factor(df$filter, rev(df$filter))

ggplot(df)+
  geom_bar(aes(filter,-sample_num_remaining), stat="identity", fill="grey70", color="black")+
  geom_bar(aes(filter,sample_num_lost), stat="identity",fill="darkred", color="black")+
  geom_text(aes(x=filter, y=-min(sample_num_remaining)/2,  label=comma(sample_num_remaining)))+
  geom_text(aes(x=filter, y=max(sample_num_lost)/1.5,  label=comma(sample_num_lost)))+
  geom_hline(yintercept=0)+
  coord_flip()+theme_bw()+ylab("")+xlab("")+
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "grey20", size=12),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  scale_x_discrete(position = "top")

ggsave(here("DNAm/figs","thirdbatch_probe_attrition.pdf"), width = 8, height = 3)
ggsave(here("DNAm/figs/jpeg","thirdbatch_probe_attrition.jpeg"), width = 8, height = 3)

epic.organoid<-epic.organoid[,which(!(colnames(epic.organoid)%in%c("array.id.path","plate_path","GEPadGI","BOX","Source","Sample_Group", "Pool_ID")))]

thirdbatch_organoid_beta<-organoid_beta
thirdbatch_epic.organoid<-epic.organoid

save(thirdbatch_organoid_beta, thirdbatch_epic.organoid, file=here("DNAm/data/batch_three_raw","thirdbatch_betas_normalized.RData"))








load(file=here("DNAm/data/batch_three_raw","thirdbatch_betas_normalized.RData"))
thirdbatch_epic.organoid$Segment<-as.factor(thirdbatch_epic.organoid$Segment)
levels(thirdbatch_epic.organoid$Segment)<-c("NA","DUO","SC","TI")

## tidier AZA
thirdbatch_epic.organoid$AZA_patient_treated<-thirdbatch_epic.organoid$AZA



## passage numeric
thirdbatch_epic.organoid$passage.or.rescope.no_numeric<-as.factor(thirdbatch_epic.organoid$passage.or.rescope.no)
levels(thirdbatch_epic.organoid$passage.or.rescope.no_numeric)<-c(1,10,12,2,3,4,6,7,8,9)
thirdbatch_epic.organoid$passage.or.rescope.no_numeric<-as.numeric(as.character(thirdbatch_epic.organoid$passage.or.rescope.no_numeric))





#' ### Principal Component Analysis (PCA)
pca_res <- prcomp(t(thirdbatch_organoid_beta))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

thirdbatch_epic.organoid$Sentrix_ID<-as.factor(as.character(thirdbatch_epic.organoid$sentrix.id))

meta_categorical <- thirdbatch_epic.organoid[, c(5,8,9,11,16,19)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(thirdbatch_epic.organoid[, c(10, 45)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Sentrix ID","Diagnosis","Gender","Segment","Age Group","Treatment")
colnames(meta_continuous) <- c( "age","passage_numeric")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8))

ggsave(here("DNAm/figs/thirdbatch_heat_scree.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)
ggsave(here("DNAm/figs/jpeg","thirdbatch_heat_scree.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)




## PC vs PC plot
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, thirdbatch_epic.organoid, by="array.id")
Loadings_meta$Segment<-as.character(Loadings_meta$Segment)

ggplot(Loadings_meta, aes(PC1, PC2, fill=Segment))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (30%)")+ylab("PC2 (20%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_sampsite
ggsave(here("DNAm/figs","thirdbatch_PC12_site.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","thirdbatch_PC12_site.jpeg"), width = 5, height = 4)



ggplot(Loadings_meta, aes(PC1, PC2, fill=age.group))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (30%)")+ylab("PC2 (20%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))
ggsave(here("DNAm/figs","thirdbatch_PC12_age.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","thirdbatch_PC12_age.jpeg"), width = 5, height = 4)


ggplot(Loadings_meta, aes(PC1, PC2, fill=treatment))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (30%)")+ylab("PC2 (20%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))
ggsave(here("DNAm/figs","thirdbatch_PC12_treatment.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","thirdbatch_PC12_treatment.jpeg"), width = 5, height = 4)







ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage.or.rescope.no_numeric)))+
  geom_point(shape=21,size=3, color="black")+
  scale_fill_manual(values=pass_col, name="Passage")+theme_bw()+
  xlab("PC2 (20%)")+ylab("PC3 (10%)")+th+theme(axis.text = element_text(size=12),
                                                axis.title = element_text(size=14),
                                                plot.margin = margin(0.6, 1, 0.6, 1, "cm"))

ggsave(here("DNAm/figs","thirdbatch_PC23_passage_confoundedthough.pdf"), width = 5.5, height = 4)
ggsave(here("DNAm/figs/jpeg","thirdbatch_PC23_passage_confoundedthough.jpeg"), width = 5.5, height = 4)



ggplot(Loadings_meta, aes(PC1, PC2, fill=diagnosis))+
  geom_point(shape=21,size=3, color="black")+
  fillscale_diagnosis+
  theme_bw()+
  xlab("PC1 (30%)")+ylab("PC2 (20%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 1, 1, 1, "cm"))
ggsave(here("DNAm/figs","thirdbatch_PC34_diagnosis.pdf"), width = 6, height = 4)
ggsave(here("DNAm/figs/jpeg","thirdbatch_PC34_diagnosis.jpeg"), width = 6, height = 4)




#'## R Session Info
sessionInfo()
