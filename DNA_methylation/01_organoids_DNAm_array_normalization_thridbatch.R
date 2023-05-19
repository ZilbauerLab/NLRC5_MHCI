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

#######################################################################
## FIRST BATCH
#######################################################################
path<-"../data/raw/all_methylation_share/ALL_IDATS"


# cohort data Fliss shared
load(here("../data/raw/all_methylation_share", "AllMethylationArraySamples_31Jul19_Phenotypes.RData"))


#'### Organoids are all on EPIC
#' I will use the EPIC annotation file downloaded at:
#' ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip

epic.organoid$sample_ID<-sapply(1:nrow(epic.organoid), function(x) paste(epic.organoid$case.no[x],"_",epic.organoid$sample.site[x], sep=""))
tapply(epic.organoid$case.no, list(epic.organoid$diagnosis), function(x) length(unique(x)))
tapply(epic.organoid$sample_ID, list(epic.organoid$sample.site, epic.organoid$diagnosis), function(x) length(unique(x)))
duplicates<-tapply(epic.organoid$sample_ID, list(epic.organoid$sample.site, epic.organoid$diagnosis), length)-tapply(epic.organoid$sample_ID, list(epic.organoid$sample.site, epic.organoid$diagnosis), function(x) length(unique(x)))
#Samples are duplicated as multiple passages of the same line are included
tapply(epic.organoid$sample_ID, list(epic.organoid$sample.site, epic.organoid$passage.or.rescope.no), function(x) length(unique(x)))

repeated_passage<-names(table(epic.organoid$sample_ID)[which(table(epic.organoid$sample_ID)>1)])
repeated_passage<-epic.organoid[which(epic.organoid$sample_ID%in%repeated_passage),]
tapply(repeated_passage$sample_ID, list(repeated_passage$case.no,repeated_passage$passage.or.rescope.no, repeated_passage$sample.site), function(x) length(x))





#' # Preprocess
epic.organoid$sample_ID<-sapply(1:nrow(epic.organoid), function(x) paste(epic.organoid$case.no[x],"_",epic.organoid$sample.site[x], sep=""))
epic.organoid$sentrix_ID<-gsub("\\_.*","",epic.organoid$array.id)
epic.organoid$array.id.path <- file.path(here(path, epic.organoid$array.id))

rgset_ibd_EPIC_organoid <- read.metharray(epic.organoid$array.id.path, verbose = TRUE,force=TRUE)# multiple DMAP files common with epic so need to force https://support.bioconductor.org/p/97773/

#' Background and dye bias correction with noob thhrough funnorm implemented in minfi
#' http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_ibd_EPIC_organoid, sex=epic.organoid$sex)

ibd_EPIC_organoid_beta<-getBeta(MSet.illumina)
dim(ibd_EPIC_organoid_beta)# n=28


#' ### Detection pvalue analysis
#' Same array as in primary has bad detection pvalue so will need to be removed, but it is only one organoid sample. sentrix 201172560041
avg_detPval <- colMeans(detectionP(rgset_ibd_EPIC_organoid))

epic.organoid$det_pval<-avg_detPval

ggplot(epic.organoid)+geom_boxplot(aes(as.factor(sentrix_ID), det_pval, fill=as.factor(sentrix_ID)), outlier.shape = NA)+
  geom_point(aes(as.factor(sentrix_ID), det_pval, group=sample_ID, fill=as.factor(sentrix_ID)), shape=21, color="black",
             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE)

ggsave(here("DNAm/figs","detection_pvalue_organoids.pdf"), width=6, height=5)
ggsave(here("DNAm/figs/jpeg","detection_pvalue_organoids.jpeg"), width=6, height=5)


#' remove bad array until it can be sorted (one whole array, one sample)
epic.organoid<-epic.organoid[which(epic.organoid$det_pval<0.005),]


#' normalize raw again but this time without bad samples
rgset_ibd_EPIC_organoid <- read.metharray(epic.organoid$array.id.path, verbose = TRUE,force=TRUE)# multiple DMAP files common with epic so need to force https://support.bioconductor.org/p/97773/
MSet.illumina <- preprocessFunnorm(rgset_ibd_EPIC_organoid, sex=epic.organoid$sex)
ibd_EPIC_organoid_beta<-getBeta(MSet.illumina)


#' ###  beta distribution before and after normalization
set.seed(1)

#extract beta values for plotting
ibd_beta_raw<-getBeta(rgset_ibd_EPIC_organoid)
#remove 6 from raw for plot
ibd_beta_raw<-ibd_beta_raw[,which(colnames(ibd_beta_raw)%in%epic.organoid$array.id)]
identical(colnames(ibd_beta_raw),epic.organoid$array.id)

## for speed with only plot a sample of probes
Beta_melted<- melt(ibd_EPIC_organoid_beta[sample(nrow(ibd_EPIC_organoid_beta),(nrow(ibd_EPIC_organoid_beta)/subset)),])
Beta_melted_raw<- melt(ibd_beta_raw[sample(nrow(ibd_beta_raw),(nrow(ibd_beta_raw)/subset)),])

#remove NAs before plotting (otherwise get many non-inifnite warnings)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
Beta_Plot_raw<-Beta_melted_raw[which(!(is.na(Beta_melted_raw$value))),]

#add meta
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,epic.organoid, by.x="ID", by.y="array.id")
colnames(Beta_Plot_raw)<-c("CpG","ID","Beta")
Beta_Plot_raw<-merge(Beta_Plot_raw,epic.organoid, by.x="ID", by.y="array.id")


ggplot(Beta_Plot, aes(Beta, group=as.character(pdata.id), color=as.character(sample.site)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")

ggsave(here("DNAm/figs","beta_distribution_EPIC_organoid.pdf"), w=10, h=5)
ggsave(here("DNAm/figs/jpeg","beta_distribution_EPIC_organoid.jpeg"), w=10, h=5)

ggplot(Beta_Plot_raw, aes(Beta, group=as.character(pdata.id), color=as.character(sample.site)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")

ggsave(here("DNAm/figs","beta_distribution_raw_EPIC_organoid.pdf"), w=10, h=5)
ggsave(here("DNAm/figs/jpeg","beta_distribution_raw_EPIC_organoid.jpeg"), w=10, h=5)






#' ### Clustering By Any Meta Data variable 

# remove rows with NAs
ibdBetas_cluster<-ibd_EPIC_organoid_beta[complete.cases(ibd_EPIC_organoid_beta),]

d <- dist(t(ibdBetas_cluster))
hc <- hclust(d, method = "complete") #single, complete, average, ward

pdf(here("DNAm/figs","cluster_wholeEPIC_organoid.pdf"),  width=30)
myplclust(hc, labels=epic.organoid$sample_ID, lab.col=as.fumeric(epic.organoid$sample.site), cex=1.5)
dev.off()

#' #### Genotyping Probes
SNPs <- getSnpBeta(rgset_ibd_EPIC_organoid)
SNPs<-SNPs[complete.cases(SNPs),]# 65 cause one was all NA

SNPs<-SNPs[,which(colnames(SNPs)%in%epic.organoid$array.id)]
identical(colnames(SNPs),epic.organoid$array.id)

d <- dist(t(SNPs))
hc <- hclust(d, method = "complete") #single, complete, average, ward

pdf(here("DNAm/figs","cluster_snps_EPIC_organoid.pdf"),  width=30)
myplclust(hc, labels=epic.organoid$sample_ID, lab.col=as.fumeric(as.character(epic.organoid$case.no)), cex=1.5)
dev.off()

#' #### Sex clustering
#' Using the cg ID to chromosome annotation from illumina 
#' https://emea.support.illumina.com/downloads/infinium-methylationepic-v1-0-product-files.html
anno_EPIC<-read.csv(here("../data", "MethylationEPIC_v-1-0_B4.csv"), skip=7)
ibd_EPIC_organoid_beta<-ibd_EPIC_organoid_beta[which(rownames(ibd_EPIC_organoid_beta)%in%anno_EPIC$IlmnID),]

#' save with sex chromosome for seperate analysis
save(ibd_EPIC_organoid_beta, file=paste(here("../data"),"/ibd_EPIC_organoid_beta_sex.RData",sep=""))

anno_EPIC<-anno_EPIC[match(rownames(ibd_EPIC_organoid_beta),anno_EPIC$IlmnID),]
identical(rownames(ibd_EPIC_organoid_beta),anno_EPIC$IlmnID)

ibd_sex<-ibd_EPIC_organoid_beta[which(anno_EPIC$CHR%in%c('X','Y')),]

d <- dist(t(ibd_sex))
hc <- hclust(d, method = "complete") #single, complete, average, ward

pdf(here("DNAm/figs","cluster_sex_EPIC_organoid.pdf"),  width=30)
myplclust(hc, labels=epic.organoid$sample_ID, lab.col=as.fumeric(epic.organoid$sex), cex=1.5)
dev.off()

pdf(here("DNAm/figs","cluster_sex_EPIC_organoid_arrayID.pdf"),  width=30)
myplclust(hc, labels=epic.organoid$array.id, lab.col=as.fumeric(epic.organoid$sex), cex=1.5)
dev.off()


epic.organoid<-epic.organoid[which(epic.organoid$array.id%in%colnames(ibd_sex)),]
identical(colnames(ibd_sex), as.character(epic.organoid$array.id))
d <- dist(t(ibd_sex))
hc <- hclust(d, method = "complete") #single, complete, average, ward
myplclust(hc, labels=epic.organoid$sample_ID, lab.col=as.fumeric(epic.organoid$sex), cex=1.5)


#' #### Remove samples which do not cluster correctly
#' T036 is  sex mismatch 287 is a sample site mixup
epic.organoid<-epic.organoid[which(!(epic.organoid$array.id%in%c("203548970031_R03C01","203548970036_R03C01"))),]
ibd_EPIC_organoid_beta<-ibd_EPIC_organoid_beta[,which(!(colnames(ibd_EPIC_organoid_beta)%in%c("203548970031_R03C01","203548970036_R03C01")))]
identical(colnames(ibd_EPIC_organoid_beta), as.character(epic.organoid$array.id))
dim(ibd_EPIC_organoid_beta)



#' ### Probe Filtering 
# SNP probes should already be removed
ibd_EPIC_organoid_beta <- ibd_EPIC_organoid_beta[!grepl("rs",rownames(ibd_EPIC_organoid_beta)), ]
dim(ibd_EPIC_organoid_beta) # probes = 865859, n = 168

#' #### Sex Chromosomes 
anno_EPIC<-anno_EPIC[anno_EPIC$IlmnID%in%rownames(ibd_EPIC_organoid_beta),]
identical(rownames(ibd_EPIC_organoid_beta),anno_EPIC$IlmnID)
ibd_EPIC_organoid_beta <- ibd_EPIC_organoid_beta[!anno_EPIC$CHR%in%c("X", "Y"), ]
dim(ibd_EPIC_organoid_beta) # probes = 846232, n = 168, 19627 filtered
filt_sex<-nrow(ibd_EPIC_organoid_beta)


#' #### Cross-hybridizing probes and polymorphic probes. 
#' https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1
#' "43,254 cross-reactive probes with ≥ 47 bp homology with an off-target site, of which 15,782 (36.5 %) are new to the EPIC platform"
#' They include this annotated list in their supplement.
#' wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM2_ESM.csv
cross_reactive<-read.csv(here("../data","13059_2016_1066_MOESM2_ESM.csv"), stringsAsFactors = F)
ibd_EPIC_organoid_beta<-ibd_EPIC_organoid_beta[which(!(rownames(ibd_EPIC_organoid_beta)%in%cross_reactive$PROBE)),]
dim(ibd_EPIC_organoid_beta) # probes = 814600, n = 168, 31632 filtered
filt_cross<-nrow(ibd_EPIC_organoid_beta)


#For polymorphic probes I will The Pidsley annotation aswell for "Probes overlapping genetic variants at targeted CpG sites." 
#and "Probes overlapping genetic variants at single base extension sites for Infinium Type I probes" 
#but NOT "Probes with genetic variants overlapping the body of the probe: 48 base pairs for Infinium Type I probes and 49 base pairs for Infinium Type II probes."

#wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM4_ESM.csv
#wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM5_ESM.csv

polymorphic<-read.csv(here("../data","13059_2016_1066_MOESM4_ESM.csv"),stringsAsFactors = F)
length(unique(polymorphic$PROBE))
baseext<-read.csv(here("../data","13059_2016_1066_MOESM5_ESM.csv"),stringsAsFactors = F)
length(unique(baseext$PROBE))
ibd_EPIC_organoid_beta<-ibd_EPIC_organoid_beta[which(!(rownames(ibd_EPIC_organoid_beta)%in%c(polymorphic$PROBE, baseext$PROBE))),]
dim(ibd_EPIC_organoid_beta) # probes = 802848, n = 168, 11752 filtered
filt_poly<-nrow(ibd_EPIC_organoid_beta)












#' ## probe filtering detection pvalue and backgrund

#' Remove probes with high NA count
na_count_probe <-sapply(1:nrow(ibd_EPIC_organoid_beta), function(y) length(which(is.na(ibd_EPIC_organoid_beta[y,]))))
na_count_probe_good<-which(na_count_probe<(ncol(ibd_EPIC_organoid_beta)*0.05))
ibd_EPIC_organoid_beta<-ibd_EPIC_organoid_beta[na_count_probe_good,]
dim(ibd_EPIC_organoid_beta)# probes = 802848, n = 168, 0 filtered
filt_bead<-nrow(ibd_EPIC_organoid_beta)


#' detection pval
detP <- detectionP(rgset_ibd_EPIC_organoid)
detP<-detP[,which(colnames(detP)%in%epic.organoid$array.id)]
identical(colnames(detP),epic.organoid$array.id)

failed <- detP>0.05
bad_det_p<-names(which(rowMeans(failed)>0.01))
bad_det_psamp<-names(which(colMeans(failed)>0.01))

ibd_EPIC_organoid_beta<-ibd_EPIC_organoid_beta[which(!(rownames(ibd_EPIC_organoid_beta)%in%bad_det_p)),]
ibd_EPIC_organoid_beta<-ibd_EPIC_organoid_beta[,which(!(colnames(ibd_EPIC_organoid_beta)%in%bad_det_psamp))]

dim(ibd_EPIC_organoid_beta)# probes = 800302, n = 168,  35409 filtered
filt_detp<-nrow(ibd_EPIC_organoid_beta)


save(ibd_EPIC_organoid_beta, epic.organoid, file=paste(here("../data"),"/ibd_beta_organoids.RData",sep="")) 


#######################################################################
## SECOND BATCH
#######################################################################
#' ### Load Data
path<-"DNAm/data/newbatch"
epic.organoid1<-read.csv(here(path, "METH_K_Nayak_20035_SS_Plate1.csv"), skip=7)
epic.organoid1$plate_path<-paste(epic.organoid1$Sample_Plate,"_96_Samples",sep="")
epic.organoid2<-read.csv(here(path, "METH_K_Nayak_20035_SS_Plate2.csv"), skip=7)
epic.organoid2$plate_path<-paste(epic.organoid2$Sample_Plate,"_96_samples",sep="")
epic.organoid3<-read.csv(here(path, "METH_K_Nayak_20035_SS_Plate3.csv"), skip=7)
epic.organoid3$plate_path<-paste(epic.organoid3$Sample_Plate,"_16_samples",sep="")

epic.organoid<-rbind(epic.organoid1,epic.organoid2,epic.organoid3)

#'### add gender
sampleInfo<-read.csv(here("DNAm/data/newbatch", "Final samples for submission 2021.csv"))
epic.organoid<-merge(epic.organoid,sampleInfo[c("Sample.Name", "Age", "Diagnosis", "Gender", "GEPadGI", "BOX", "Segment","Source", "Wnt.type","X.1","Biobank.Rachel.Replicates")], by.x="Sample_Name", by.y="Sample.Name")

#' ### Normalize DNAm Arrays
here(path)
epic.organoid$array.id.path <- file.path(here(path,"Data"),epic.organoid$plate_path, epic.organoid$Sentrix_ID, paste(epic.organoid$Sentrix_ID, epic.organoid$Sentrix_Position, sep="_"))
epic.organoid$array.id<-paste(epic.organoid$Sentrix_ID, epic.organoid$Sentrix_Position, sep="_")
epic.organoid$individual<-sapply(1:nrow(epic.organoid), function(x) strsplit(epic.organoid$Sample_Name[x]," ")[[1]][1])

epic.organoid$passage<-as.numeric(sapply(1:nrow(epic.organoid), function(x) {
  parts<-strsplit(epic.organoid$Sample_Name[x]," ")[[1]]
  parts<-parts[which(!(parts=="WFP"))]
  gsub("P|p","",parts[grep("P|p", parts)])}))
epic.organoid$passage_hilo<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$passage[x]<5){"low"}else{"high"})
epic.organoid$condition<-sapply(1:nrow(epic.organoid), function(x) strsplit(epic.organoid$Sample_Name[x]," ")[[1]][4])
epic.organoid$condition<-as.factor(epic.organoid$condition)
levels(epic.organoid$condition)<-c("D","FB","IFNg","IFNg",NA,"TNFa","UD","UT",NA)
epic.organoid$condition<-as.character(epic.organoid$condition)
epic.organoid$comparison<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$condition[x]%in%c("UD","D")){"differentiation"}else{"cytokine"})
epic.organoid$treatment<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$comparison[x]=="cytokine"){epic.organoid$condition[x]}else{"UT"})
epic.organoid$differentiation<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$comparison[x]=="differentiation"){epic.organoid$condition[x]}else{"UD"})
colnames(epic.organoid)[17]<-"diagnosis_simplified"


# multiple DMAP files common with epic so need to force https://support.bioconductor.org/p/97773/
rgset_organoid <- read.metharray(epic.organoid$array.id.path, verbose = FALSE,force=TRUE)

# Background and dye bias correction with noob thhrough funnorm implemented in minfi
#http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_organoid, sex=epic.organoid$Gender)
organoid_beta<-getBeta(MSet.illumina)

print(paste("Samples available: ",ncol(organoid_beta),"; Probes available: ",nrow(organoid_beta),sep=""))

save(rgset_organoid, MSet.illumina,epic.organoid, file=here("DNAm/data/newbatch","newbatch_raw.RData"))





load(here("DNAm/data/newbatch","newbatch_raw.RData"))
organoid_beta<-getBeta(MSet.illumina)

#' ### Detection p values across all probes for each sample
avg_detPval <- colMeans(detectionP(rgset_organoid))
epic.organoid$det_pval<-avg_detPval

ggplot(epic.organoid)+geom_boxplot(aes(as.factor(Sentrix_ID), det_pval, fill=as.factor(Sentrix_ID)), outlier.shape = NA)+
  geom_point(aes(as.factor(Sentrix_ID), det_pval, group=Sample_Name, fill=as.factor(Sentrix_ID)), shape=21, color="black",
             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE)

ggsave(here("DNAm/figs","newbatch_detection_pvalue_organoids.pdf"), width=6, height=5)
ggsave(here("DNAm/figs/jpeg","newbatch_detection_pvalue_organoids.jpeg"), width=6, height=5)




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

ggsave(here("DNAm/figs","newbatch_beta_distribution_organoid.pdf"),grid.arrange(beta_dis_EPIC_raw,beta_dis_EPIC),  w=10, h=5)
ggsave(here("DNAm/figs/jpeg","newbatch_beta_distribution_organoid.jpeg"),grid.arrange(beta_dis_EPIC_raw,beta_dis_EPIC), w=10, h=5)



#'#### Confirm individuals ID with SNPs probes and clustering by DNAm

# remove rows with NAs
Betas_cluster<-organoid_beta[complete.cases(organoid_beta),]

d <- dist(t(Betas_cluster))
hc <- hclust(d, method = "complete") #single, complete, average, ward
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(epic.organoid$Segment), cex=1.5)

pdf(here("DNAm/figs","newbatch_cluster_wholeEPIC_organoid.pdf"), width=30)
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

pdf(here("DNAm/figs","newbatch_cluster_snps_EPIC_organoid.pdf"), width=30)
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

pdf(here("DNAm/figs","newbatch_cluster_sex_EPIC_organoid.pdf"), width=30)
myplclust(hc, labels=epic.organoid$Sample_Name, lab.col=as.fumeric(epic.organoid$Gender), cex=1.5)
dev.off()

#'##' Another cehck for sex mixups
epic.organoid$sex_mix<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$individual[x] %in% c("T036","T279","T203","T193","T238","T202","T196","T192","T362")){
  epic.organoid$Sample_Name[x]}else{epic.organoid$Gender[x]})

organoid_beta_sex <- organoid_beta[anno_EPIC$CHR%in%c("X", "Y"), ]
save(organoid_beta_sex, epic.organoid, file=here("DNAm/data/newbatch","newbatch_betas_normalized_sexonly.RData"))

load(file=here("DNAm/data/newbatch","newbatch_betas_normalized_sexonly.RData"))


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

ggsave(here("DNAm/figs","newbatch_beta_distribution_EPIC_Y.pdf"),grid.arrange(all_sample,sample_split), w=10, h=15)
ggsave(here("DNAm/figs/jpeg","newbatch_beta_distribution_EPIC_Y.jpeg"),grid.arrange(all_sample,sample_split), w=10, h=15)

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

ggsave(here("DNAm/figs","newbatch_beta_distribution_EPIC_X.pdf"),grid.arrange(all_sample,sample_split), w=10, h=15)
ggsave(here("DNAm/figs/jpeg","newbatch_beta_distribution_EPIC_X.jpeg"),grid.arrange(all_sample,sample_split), w=10, h=15)




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

ggsave(here("DNAm/figs","newbatch_probe_attrition.pdf"), width = 8, height = 3)
ggsave(here("DNAm/figs/jpeg","newbatch_probe_attrition.jpeg"), width = 8, height = 3)

epic.organoid<-epic.organoid[,which(!(colnames(epic.organoid)%in%c("array.id.path","plate_path","GEPadGI","BOX","Source","Sample_Group", "Pool_ID")))]

newbatch_organoid_beta<-organoid_beta
newbatch_epic.organoid<-epic.organoid

save(newbatch_organoid_beta, newbatch_epic.organoid, file=here("DNAm/data/newbatch","newbatch_betas_normalized.RData"))








load(file=here("DNAm/data/newbatch","newbatch_betas_normalized.RData"))
newbatch_epic.organoid$Segment<-as.factor(newbatch_epic.organoid$Segment)
levels(newbatch_epic.organoid$Segment)<-c("Duo","Duo","SC","TI")

#' ### Principal Component Analysis (PCA)
pca_res <- prcomp(t(newbatch_organoid_beta))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

newbatch_epic.organoid$Sentrix_ID<-as.factor(as.character(newbatch_epic.organoid$Sentrix_ID))

meta_categorical <- newbatch_epic.organoid[, c(3,4,11, 8,9,10,14,17)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(newbatch_epic.organoid[, c(6, 15)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Plate","Sentrix ID","Diagnosis","Gender","Segment","WNT Type","Individual", "passage_hilo")
colnames(meta_continuous) <- c( "age","passage_numeric")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8))

ggsave(here("DNAm/figs/newbatch_heat_scree.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)
ggsave(here("DNAm/figs/jpeg","newbatch_heat_scree.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)




## PC vs PC plot
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, newbatch_epic.organoid, by="array.id")

ggplot(Loadings_meta, aes(PC1, PC2, fill=Segment))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (33%)")+ylab("PC2 (13%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_sampsite
ggsave(here("DNAm/figs","newbatch_PC12_site.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","newbatch_PC12_site.jpeg"), width = 5, height = 4)

Loadings_meta$possible_mixed<-""
Loadings_meta$possible_mixed[which(Loadings_meta$Sample_Name%in%c("T419 DUO P2","T419 SC P1","T363 DUO P2","T311 DUO P10","T346 TI P2"))]<-Loadings_meta$Sample_Name[which(Loadings_meta$Sample_Name%in%c("T419 DUO P2","T419 SC P1","T363 DUO P2","T311 DUO P10","T346 TI P2"))]

ggplot(Loadings_meta, aes(PC1, PC2, fill=Segment))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  geom_text(aes(label=possible_mixed, color=Segment))+
  xlab("PC1 (33%)")+ylab("PC2 (13%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+fillscale_sampsite+colscale_sampsite
ggsave(here("DNAm/figs","newbatch_PC12_site_label_mix.pdf"), width = 11, height = 10)
ggsave(here("DNAm/figs/jpeg","newbatch_PC12_site_label_mix.jpeg"), width = 11, height = 10)





ggplot(Loadings_meta, aes(PC2, PC3, fill=Segment))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC2 (13%)")+ylab("PC3 (7.5%)")+th+theme(axis.text = element_text(size=12),
                                                axis.title = element_text(size=14),
                                                plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_sampsite
ggsave(here("DNAm/figs","newbatch_PC23_site.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","newbatch_PC23_site.jpeg"), width = 5, height = 4)



ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage)))+
  geom_point(shape=21,size=3, color="black")+
  scale_fill_manual(values=pass_col, name="Passage")+theme_bw()+
  xlab("PC2 (13%)")+ylab("PC3 (7.5%)")+th+theme(axis.text = element_text(size=12),
                                                axis.title = element_text(size=14),
                                                plot.margin = margin(0.6, 1, 0.6, 1, "cm"))

ggsave(here("DNAm/figs","newbatch_PC23_passage.pdf"), width = 5.5, height = 4)
ggsave(here("DNAm/figs/jpeg","newbatch_PC23_passage.jpeg"), width = 5.5, height = 4)

cor(Loadings_meta$PC3, Loadings_meta$passage, method="spearman")

Loadings_meta$diagnosis_simplified<-as.factor(Loadings_meta$diagnosis_simplified)
levels(Loadings_meta$diagnosis_simplified)<-c("CD","Control","IBD","UC")

ggplot(Loadings_meta, aes(PC3, PC4, fill=diagnosis_simplified))+
  geom_point(shape=21,size=3, color="black")+
  fillscale_diagnosis+
  theme_bw()+
  xlab("PC3 (7.5%)")+ylab("PC4 (3%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 1, 1, 1, "cm"))
ggsave(here("DNAm/figs","newbatch_PC34_diagnosis.pdf"), width = 6, height = 4)
ggsave(here("DNAm/figs/jpeg","newbatch_PC34_diagnosis.jpeg"), width = 6, height = 4)





#######################
## Relabel and removed samples
#######################
# gender mis labelled
newbatch_epic.organoid$Gender[grep("T279|T193",newbatch_epic.organoid$Sample_Name)]<-"F"
newbatch_epic.organoid$Gender[grep("T196|T192|T362",newbatch_epic.organoid$Sample_Name)]<-"M"


# Fliss note
# T202/T203: It has become clear that the sample labels for T202 & T203 were switched. The "T203" sample 
# (actually T202) was removed at initial sample QC, so we'll leave it out. the "T202" sample (actually T203) 
# was mislabelled as male, so was retained & wasn't picked up in sample QC (as T203 actually is male). 
# Relabel the case.no for this sample & update the age. Both T202 & T203 have CD
newbatch_epic.organoid[grep("202|203",newbatch_epic.organoid$Sample_Name),]

T202<-newbatch_epic.organoid[grep("202",newbatch_epic.organoid$Sample_Name),]
T203<-newbatch_epic.organoid[grep("203",newbatch_epic.organoid$Sample_Name),]

newbatch_epic.organoid[grep("203",newbatch_epic.organoid$Sample_Name),c("Sample_Well", "Sample_Plate","Sentrix_ID", "Sentrix_Position","array.id" )]<-T202[,c("Sample_Well", "Sample_Plate","Sentrix_ID", "Sentrix_Position","array.id" )]
newbatch_epic.organoid[grep("202",newbatch_epic.organoid$Sample_Name),c("Sample_Well", "Sample_Plate","Sentrix_ID", "Sentrix_Position","array.id" )]<-T203[,c("Sample_Well", "Sample_Plate","Sentrix_ID", "Sentrix_Position","array.id" )]



# 419 sample site switch
newbatch_epic.organoid$Sample_Name[grep("205605880094_R06C01",newbatch_epic.organoid$array.id)]<-"T419 SC P1"
newbatch_epic.organoid$Segment[grep("205605880094_R06C01",newbatch_epic.organoid$array.id)]<-"SC"
newbatch_epic.organoid$Sample_Name[grep("205605880096_R08C01",newbatch_epic.organoid$array.id)]<-"T419 DUO P2"
newbatch_epic.organoid$Segment[grep("205605880096_R08C01",newbatch_epic.organoid$array.id)]<-"Duo"

# age 10 months to a numeric
newbatch_epic.organoid$Age[grep("10m", newbatch_epic.organoid$Age)]<-"0.83"

# remove samples that mislabelling can't be explained
newbatch_epic.organoid<-newbatch_epic.organoid[which(!(newbatch_epic.organoid$Sample_Name%in%c("T346 TI P2","T363 DUO P2","T036 SC P11","T238 TI P6"))),]

table(newbatch_epic.organoid$Biobank.Rachel.Replicates)
newbatch_epic.organoid$sex_mix<-NULL

newbatch_organoid_beta<-newbatch_organoid_beta[,which(colnames(newbatch_organoid_beta)%in%newbatch_epic.organoid$array.id)]
newbatch_epic.organoid<-newbatch_epic.organoid[match(colnames(newbatch_organoid_beta), newbatch_epic.organoid$array.id),]
identical(colnames(newbatch_organoid_beta),newbatch_epic.organoid$array.id )

#' ### Principal Component Analysis (PCA)
pca_res <- prcomp(t(newbatch_organoid_beta))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

newbatch_epic.organoid$Sentrix_ID<-as.factor(as.character(newbatch_epic.organoid$Sentrix_ID))


meta_categorical <- newbatch_epic.organoid[, c(3,4,11, 8,9,10,14,17)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(newbatch_epic.organoid[, c(6, 15)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Plate","Sentrix ID","Diagnosis","Gender","Segment","WNT Type","Individual", "passage_hilo")
colnames(meta_continuous) <- c( "age","passage_numeric")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8))

ggsave(here("DNAm/figs/newbatch_heat_scree_aftersampleQC.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)
ggsave(here("DNAm/figs/jpeg","newbatch_heat_scree_aftersampleQC.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)




## PC vs PC plot
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, newbatch_epic.organoid, by="array.id")

ggplot(Loadings_meta, aes(PC1, PC2, fill=Segment))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (33%)")+ylab("PC2 (13%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_sampsite
ggsave(here("DNAm/figs","newbatch_PC12_site_aftersampleQC.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","newbatch_PC12_site_aftersampleQC.jpeg"), width = 5, height = 4)



ggplot(Loadings_meta, aes(PC2, PC3, fill=Segment))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC2 (13%)")+ylab("PC3 (7.5%)")+th+theme(axis.text = element_text(size=12),
                                                axis.title = element_text(size=14),
                                                plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_sampsite
ggsave(here("DNAm/figs","newbatch_PC23_site_aftersampleQC.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","newbatch_PC23_site_aftersampleQC.jpeg"), width = 5, height = 4)



ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage)))+
  geom_point(shape=21,size=3, color="black")+
  scale_fill_manual(values=pass_col, name="Passage")+theme_bw()+
  xlab("PC2 (13%)")+ylab("PC3 (7.5%)")+th+theme(axis.text = element_text(size=12),
                                                axis.title = element_text(size=14),
                                                plot.margin = margin(0.6, 1, 0.6, 1, "cm"))

ggsave(here("DNAm/figs","newbatch_PC23_passage_aftersampleQC.pdf"), width = 5.5, height = 4)
ggsave(here("DNAm/figs/jpeg","newbatch_PC23_passage_aftersampleQC.jpeg"), width = 5.5, height = 4)

cor(Loadings_meta$PC3, Loadings_meta$passage, method="spearman")

Loadings_meta$diagnosis_simplified<-as.factor(Loadings_meta$diagnosis_simplified)
levels(Loadings_meta$diagnosis_simplified)<-c("CD","Control","IBD","UC")

ggplot(Loadings_meta, aes(PC3, PC4, fill=diagnosis_simplified))+
  geom_point(shape=21,size=3, color="black")+
  fillscale_diagnosis+
  theme_bw()+
  xlab("PC3 (7.5%)")+ylab("PC4 (3%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 1, 1, 1, "cm"))
ggsave(here("DNAm/figs","newbatch_PC34_diagnosis_aftersampleQC.pdf"), width = 6, height = 4)
ggsave(here("DNAm/figs/jpeg","newbatch_PC34_diagnosis_aftersampleQC.jpeg"), width = 6, height = 4)


save(newbatch_organoid_beta, newbatch_epic.organoid, file=here("DNAm/data/newbatch","newbatch_betas_normalized_aftersampleQC.RData"))








#######################################################################
## THIRD BATCH
#######################################################################




#' ### Load Data
path<-"DNAm/data/batch_three_raw"
epic.organoid<-read.csv(here(path, "METH_K_Nayak_20038_SAMPLESHEET.csv"), skip=7)
epic.organoid$array.id<-paste(epic.organoid$Sentrix_ID, epic.organoid$Sentrix_Position, sep="_")


#'### add gender
load(here("DNAm/data", "AllEPICOrganoid_FullMetadata_23Sep22.RData"))
sampleInfo<-epic.samples[which(epic.samples$array.id%in%epic.organoid$array.id),]
epic.organoid<-merge(epic.organoid,sampleInfo, by="array.id")

#' ### Normalize DNAm Arrays
here(path)
epic.organoid$array.id.path <- file.path(here(path,"RAW_iDAT_FILES_20038_17JUN22",epic.organoid$array.id))

epic.organoid$individual<-sapply(1:nrow(epic.organoid), function(x) strsplit(epic.organoid$Sample_Name[x],"_")[[1]][1])
epic.organoid$condition<-as.factor(epic.organoid$treatment)
epic.organoid$condition<-as.character(epic.organoid$condition)
epic.organoid$comparison<-sapply(1:nrow(epic.organoid), function(x) if(epic.organoid$condition[x]%in%c("IFNg","LPS","AZA")){"UT"}else{"NT"})
epic.organoid$Segment<-epic.organoid$tissue
epic.organoid$Gender<-epic.organoid$sex


# multiple DMAP files common with epic so need to force https://support.bioconductor.org/p/97773/
rgset_organoid <- read.metharray(epic.organoid$array.id.path, verbose = FALSE,force=TRUE)

# Background and dye bias correction with noob thhrough funnorm implemented in minfi
#http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_organoid, sex=epic.organoid$Gender)
organoid_beta<-getBeta(MSet.illumina)

print(paste("Samples available: ",ncol(organoid_beta),"; Probes available: ",nrow(organoid_beta),sep=""))

save(rgset_organoid, MSet.illumina,epic.organoid, file=here("DNAm/data/batch_three_raw","thridbatch_raw.RData"))





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
