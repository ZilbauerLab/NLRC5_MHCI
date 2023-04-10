#'---
#'title: DNAm preprocessing and normalization
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' # 450K
#' ### Load Libraries
suppressMessages(library(here))
suppressMessages(library(minfi))
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(rafalib)
library(scales)

options(stringsAsFactors = FALSE)
options(scipen = 999)




#' ### Load Functions
source(here("general_functions","00_pretty_plots.R"))
source(here("general_functions","00_Heat_scree_plot_generic.R"))


#' Will be using the annotation files downloaded at
#' wget ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/ProductFiles/HumanMethylation450/HumanMethylation450_15017482_v1-2.csv
#' https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL16304 #the download table buttonMethylation Arrays



#' for quicker plotting
subset=4

#' ### Load Raw Data
#' Raw IDAT files were processed with Illumina’s GenomeStudio software and background normalised using negative control probes in minfi.
path<-"../data/raw/K450_METHYLATION_DATA"

sampleinfo <- read.table(here(path, "Zerbino_K450MethylationArray_SampleInfo.txt"), header=T, sep="\t")

sampleinfo$sample_ID<-sapply(1:nrow(sampleinfo), function(x) paste(sampleinfo$case.no[x],"_",sampleinfo$sample.site[x], sep=""))
sampleinfo$sentrix_ID<-gsub("\\_.*","",sampleinfo$array.id)

sampleinfo$array.id.path <- file.path(here(path, sampleinfo$array.id))
rgset_ibd <- read.metharray(sampleinfo$array.id.path, verbose = TRUE)



#' Background and dye bias correction with noob thhrough funnorm implemented in minfi
#' http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_ibd, sex=sampleinfo$sex)

ibd_beta<-getBeta(MSet.illumina)





#' ### Detection pvalue analysis
#' Using the control probes to estimate background on the array, an p value is generated to define the confidence in a call at a probe. Generally across 11 arrays none seem to have wide spread failure. The detection p values will be interrogated later for individual probes and sample issues. 
avg_detPval <- colMeans(detectionP(rgset_ibd))

sampleinfo$det_pval<-avg_detPval

ggplot(sampleinfo)+geom_boxplot(aes(as.factor(sentrix_ID), det_pval, fill=as.factor(sentrix_ID)), outlier.shape = NA)+
  geom_point(aes(as.factor(sentrix_ID), det_pval, group=sample_ID, fill=as.factor(sentrix_ID)), shape=21, color="black",
             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE)+ylim(0,0.008)

ggsave(here("DNAm/figs","detection_pvalue_450k.pdf"), width=6, height=5)
ggsave(here("DNAm/figs/jpeg","detection_pvalue_450k.jpeg"), width=6, height=5)


# remove any arrays with systematically bad detection pvalues (none in 450k though)
sampleinfo<-sampleinfo[which(sampleinfo$det_pval<0.005),]
ibd_beta<-ibd_beta[,which(colnames(ibd_beta)%in%sampleinfo$array.id)]
identical(colnames(ibd_beta),sampleinfo$array.id)





#' ###  beta distribution before and after normalization
#' Samples do not seem to have a drastically different beta distribution based on sample site. Will use quantro to confirm statistically before any between quantile normalisation. On the left is the raw beta distributions and the left is after functional normalization.
set.seed(1)
ibd_beta_raw<-getBeta(rgset_ibd)

#extract beta values for plotting
ibd_beta<-getBeta(MSet.illumina)

## for speed with only plot a sample of probes
Beta_melted<- melt(ibd_beta[sample(nrow(ibd_beta),(nrow(ibd_beta)/subset)),])
Beta_melted_raw<- melt(ibd_beta_raw[sample(nrow(ibd_beta_raw),(nrow(ibd_beta_raw)/subset)),])

#remove NAs before plotting (otherwise get many non-inifnite warnings)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
Beta_Plot_raw<-Beta_melted_raw[which(!(is.na(Beta_melted_raw$value))),]

#add meta
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,sampleinfo, by.x="ID", by.y="array.id")
colnames(Beta_Plot_raw)<-c("CpG","ID","Beta")
Beta_Plot_raw<-merge(Beta_Plot_raw,sampleinfo, by.x="ID", by.y="array.id")


ggplot(Beta_Plot, aes(Beta, group=as.character(pdata.id), color=as.character(sample.site)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")

ggsave(here("DNAm/figs","beta_distribution_450K.pdf"), w=10, h=5)
ggsave(here("DNAm/figs/jpeg","beta_distribution_450K.jpeg"), w=10, h=5)

ggplot(Beta_Plot_raw, aes(Beta, group=as.character(pdata.id), color=as.character(sample.site)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")

ggsave(here("DNAm/figs","beta_distribution_raw_450K.pdf"), w=10, h=5)
ggsave(here("DNAm/figs/jpeg","beta_distribution_raw_450K.jpeg"), w=10, h=5)





#' ### Clustering By Any Meta Data variable 
#' Based on clustering of the entire 450K samples group by sample site. There is a homogenous TI cluster, with some TI grouping elsewhere, AC and SC form heterogeneous mixes. One sample (44_SC) groups away from the others and is a potential outlier.

# remove rows with NAs
ibdBetas_cluster<-ibd_beta[complete.cases(ibd_beta),]

d <- dist(t(ibdBetas_cluster))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(sampleinfo$sample.site), cex=1.5)

pdf(here("DNAm/figs","cluster_whole450k.pdf"), width=30)
myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(sampleinfo$sample.site), cex=1.5)
dev.off()

hclusters <- cutree(hc, h=78)
table(true=sampleinfo$sample.site, cluster=hclusters)



#' Genotyping Probes
#' Confirm ID with SNPs
#' Based on clustering the samples only on the 64 SNP probes on the 450k, samples should pair by patient. As seen below samples cluster perfectly by patient.

SNPs <- getSnpBeta(rgset_ibd)
SNPs<-SNPs[complete.cases(SNPs),]# 65 cause one was all NA
SNPs_450K<-SNPs[complete.cases(SNPs),]# for save later to match to genotype data

d <- dist(t(SNPs))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(as.character(sampleinfo$case.no)), cex=1.5)

pdf(here("DNAm/figs","cluster_snps.pdf"), width=30)
myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(as.character(sampleinfo$case.no)), cex=1.5)
dev.off()


#' by sex too
#' Based on clustering the samples on the X and Y chromosomes, samples cluster perfectly by reported sex. 
anno_450k<-read.csv(here("DNAm/data","HumanMethylation450_15017482_v1-2.csv"), skip=7)
anno_450k<-anno_450k[match(rownames(ibd_beta),anno_450k$IlmnID),]

ibd_sex<-ibd_beta[which(anno_450k$CHR%in%c('X','Y')),]
d <- dist(t(ibd_sex))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(sampleinfo$sex), cex=1.5)

pdf(here("DNAm/figs","cluster_sex.pdf"), width=30)
myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(sampleinfo$sex), cex=1.5)
dev.off()




#' #### probe filtering
#' SNP probes already filtered by overlapping with 450K annotation. Sex chromosomes have a very different distribution between sexes so they will be removed prior to normalization. Probes shown to hybridze to multiple parts of the genome will be filtered. Probes with a known SNP at the C of G of the queried CpG will be filtered. Probes with a high number of NA resulting from low beadcount (in 5% of samples) will be removed. Probes with a high detection p value (>0.05) in 1% of samples will be removed. The number of filtered probes at each stage is shown in the plot below.

#' ## Sex CHR
anno_450k<-anno_450k[match(rownames(ibd_beta),anno_450k$IlmnID),]

ibd_beta<-ibd_beta[which(!(anno_450k$CHR%in%c('X','Y'))),] #485512
dim(ibd_beta) # probes = 473864, n = 115, 11648 filtered
filt_sex<-nrow(ibd_beta)

#' ##Cross hybridizing probes
#' Some probes have been found to cross-hybridize with other chromosomes (Price et al. 2013 Epigenetics). Removing here
price<-read.table(here("../data","GPL16304-47833.txt"), sep='\t', header=T, skip=22)
price<-price[match(rownames(ibd_beta),price$ID),]

cross_hyb<-price[which(price$XY_Hits=="XY_YES" | price$Autosomal_Hits=="A_YES"),]
ibd_beta<-ibd_beta[which(!(rownames(ibd_beta)%in%cross_hyb$ID)),]
dim(ibd_beta) # probes = 433274, n = 115, 40590 filtered
filt_cross<-nrow(ibd_beta)


#' ## Polymorphic probes
SnpatCpG<-price[which(price$Target.CpG.SNP!=""),] # 20696
ibd_beta<-ibd_beta[which(!(rownames(ibd_beta)%in%SnpatCpG$ID)),]
dim(ibd_beta) # probes = 415080, n = 115, 18194 filtered
filt_poly<-nrow(ibd_beta)



#' ## probe filtering detection pvalue and backgrund

#' Remove probes with high NA count
na_count_probe <-sapply(1:nrow(ibd_beta), function(y) length(which(is.na(ibd_beta[y,]))))
na_count_probe_good<-which(na_count_probe<(ncol(ibd_beta)*0.05))
ibd_beta<-ibd_beta[na_count_probe_good,]
dim(ibd_beta)# probes = 415080, n = 115, 0 filtered
filt_bead<-nrow(ibd_beta)


#' detection pval
detP <- detectionP(rgset_ibd)
failed <- detP>0.05
bad_det_p<-names(which(rowMeans(failed)>0.01))
bad_det_psamp<-names(which(colMeans(failed)>0.01))

ibd_beta<-ibd_beta[which(!(rownames(ibd_beta)%in%bad_det_p)),]
ibd_beta<-ibd_beta[,which(!(colnames(ibd_beta)%in%bad_det_psamp))]

dim(ibd_beta)# probes = 414293, n = 115, 787 filtered
filt_detp<-nrow(ibd_beta)


#' these numbers match pfilter when run on raw rgset data
#' 0 samples having 1 % of sites with a detection p-value greater than 0.05 were removed
#' 18 sites were removed as beadcount <3 in 5 % of samples
#' 769 sites having 1 % of samples with a detection p-value greater than 0.05 were removed

ibd_beta_450K<-ibd_beta

#' ### Probe attrition plot
df<-data.frame(sample_num_remaining=c(485577,485512,filt_sex,filt_cross,filt_poly,filt_bead,filt_detp),
               filter=c("450K Probe Number","Removal of SNP Probes","Removal of X and Y chromosome probes",
                        "Removal of Cross Reactive Probes",
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

ggsave(here("DNAm/figs","probe_attrition_450K.pdf"), width = 8, height = 3)
ggsave(here("DNAm/figs/jpeg","probe_attrition_450K.jpeg"), width = 8, height = 3)








#' #### PCA
#' To look at overall contributors to 450k variation PCA was run and PC loading were associated to meta data variable. Generally the sample site of the epithelial cells is the primary driver of variation, followed by the case no. (ie patient). Diagnosis comes up in PC2 but so does sample site so that could be related to some unbalance between sample site and diagnosis number. Interestingly inflammation and sample site are overlapping in PC1.
# PCA
pca_res <- prcomp(t(ibd_beta))

Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)


#Restructure meta
sampleinfo$sentrix_ID<-as.factor(sampleinfo$sentrix_ID)
sampleinfo$case.no<-as.factor(sampleinfo$case.no)

meta_categorical <- sampleinfo[, c(1,4,5,10,11,18)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(sampleinfo[, c(12)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Case No.", "Diagnosis", "Sample Site","Inflammation","Sex","Sentrix ID")
colnames(meta_continuous) <- c("Age")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7))

ggsave(here("DNAm/figs","heat_scree.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7)),width = 9, height = 6)
ggsave(here("DNAm/figs/jpeg","heat_scree.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7)),width = 9, height = 6)

## PC vs PC plot
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, sampleinfo, by="array.id")

ggplot(Loadings_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=2, color="black")+theme_bw()
ggsave(here("DNAm/figs","PC1_PC2.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2.jpeg"), width = 7.5, height = 6)




######################################## meta data correlation
#' Therefore I checked the relation between meta data variables. Most (all but age) are categorical so chi-square test p values were generated and with age anova between all categorical variables.  Looks like cases were batched on arrays, inflammation associates to diagnosis (only IBD are inflammed, but not all IBD inflammed 50-60% of IBD inflammed) and sample site. Generally the sigmoid colon is more commonly inflamed (44% of SC is inflammed compared to, 29% of AC and 16%TI). "The inflammation status of a sample (inflamed vs non-inflamed) was based on the histology of a paired sample taken within 2 cm of samples at the time of the initial endoscopy."
meta_categorical <- sampleinfo[, c(1,4,5,10,11,18)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(sampleinfo[, c(12)] ) # input column numbers in meta that contain continuous variables
colnames(meta_continuous) <- c("Age")


Meta_correlation<-lapply(1:ncol(meta_continuous), function(x) as.numeric(meta_continuous[,x]))
Meta_correlation<-as.data.frame(do.call(cbind, Meta_correlation))
colnames(Meta_correlation)<-colnames(meta_continuous)

aov_meta <- sapply(1:ncol(meta_categorical), function(cat) {
  summary(aov(meta_continuous$Age ~ meta_categorical[, cat]))[[1]]$"Pr(>F)"[1]}
  )

chi_meta<-lapply(1:ncol(meta_categorical), function(cat1) {

  sapply(1:ncol(meta_categorical), function(cat2){
    suppressWarnings(chisq.test(table(meta_categorical[,cat1], meta_categorical[,cat2]))$p.value)})

  })

chi_meta<-as.data.frame(do.call(rbind, chi_meta))

colnames(chi_meta)<-colnames(meta_categorical)

#add continious
chi_meta$age<-aov_meta
chi_meta[nrow(chi_meta)+1,]<-c(aov_meta,0)

chi_meta$id<-c(colnames(meta_categorical),"age")

chi_meta<-melt(chi_meta)

chi_meta$Pvalue<-sapply(1:nrow(chi_meta), function(x) if(chi_meta$value[x]<=0.001){"<=0.001"}else{
  if(chi_meta$value[x]<=0.01){"<=0.01"}else{
    if(chi_meta$value[x]<=0.05){"<=0.05"}else{">0.05"}}})
chi_meta$variable<-as.character(chi_meta$variable)

ggplot(chi_meta, aes(id, variable, fill = Pvalue)) +
  geom_tile(color = "black",size=0.5) +  geom_text(aes(label=round(value, 3))) +
  theme_gray(8)+theme(axis.text = element_text(size =10, color="black"),
                      axis.title =  element_blank(),
                      axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5),
                      legend.text = element_text(size =10),
                      legend.title = element_text(size =10),
                      panel.background = element_blank(),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank())+
  scale_fill_manual(values=c("#084594","#4292c6","#9ecae1","#deebf7"))

ggsave(here("DNAm/figs","meta_correlation.pdf"), width = 6, height = 5)
ggsave(here("DNAm/figs/jpeg","meta_correlation.jpeg"), width = 6, height = 5)



#' ### proportion inflamed
pro_inflamed<-as.data.frame.matrix(table(sampleinfo$sample.site, sampleinfo$inflammation))
pro_inflamed$total<-rowSums(pro_inflamed)
pro_inflamed$perc_inflamed<-(pro_inflamed$Y/pro_inflamed$total)*100

pro_diagnosis<-as.data.frame.matrix(table(sampleinfo$diagnosis, sampleinfo$sample.site))
pro_diagnosis$total<-rowSums(pro_diagnosis)
pro_diagnosis$perc_AC<-(pro_diagnosis$AC/pro_diagnosis$total)*100
pro_diagnosis$perc_SC<-(pro_diagnosis$SC/pro_diagnosis$total)*100
pro_diagnosis$perc_TI<-(pro_diagnosis$TI/pro_diagnosis$total)*100













#' # EPIC
#' I will use the EPIC annotation file downloaded at:
#' ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b4-manifest-file-csv.zip

path<-"../data/raw/EPIC_METHYLATION_DATA"

sampleinfo <- read.table(here(path, "Zerbino_EPICMethylationArray_SampleInfo.txt"), header=T, sep="\t")

sampleinfo$sample_ID<-sapply(1:nrow(sampleinfo), function(x) paste(sampleinfo$case.no[x],"_",sampleinfo$sample.site[x], sep=""))
sampleinfo$sentrix_ID<-gsub("\\_.*","",sampleinfo$array.id)

sampleinfo$array.id.path <- file.path(here(path, sampleinfo$array.id))
rgset_ibd <- read.metharray(sampleinfo$array.id.path, verbose = TRUE)


#' Background and dye bias correction with noob thhrough funnorm implemented in minfi
#' http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_ibd, sex=sampleinfo$sex)
ibd_beta<-getBeta(MSet.illumina)


#'#### Detection pvalue analysis
#'Using the control probes to estimate background on the array, an p value is generated to define the confidence in a call at a probe. Across the 16 arrays one seems to have high detection p values across all samples. These samples will be removed from analysis. 
avg_detPval <- colMeans(detectionP(rgset_ibd))

sampleinfo$det_pval<-avg_detPval

ggplot(sampleinfo)+geom_boxplot(aes(as.factor(sentrix_ID), det_pval, fill=as.factor(sentrix_ID)), outlier.shape = NA)+
  geom_point(aes(as.factor(sentrix_ID), det_pval, group=sample_ID, fill=as.factor(sentrix_ID)), shape=21, color="black",
             position = position_jitter(w = 0.25))+theme_bw()+theme(axis.text.x=element_text(angle = 45, vjust = 1, hjust=1))+
  xlab("Sentrix ID")+ylab("Mean Detection P Value")+guides(fill=FALSE)

ggsave(here("DNAm/figs","detection_pvalue_EPIC.pdf"), width=6, height=5)
ggsave(here("DNAm/figs/jpeg","detection_pvalue_EPIC.jpeg"), width=6, height=5)


#' remove bad array until it can be sorted (one whole array, six samples on EPIC)
sampleinfo<-sampleinfo[which(sampleinfo$det_pval<0.005),]




#' normalize raw again but this time without bad samples
rgset_ibd <- read.metharray(sampleinfo$array.id.path, verbose = TRUE)
MSet.illumina <- preprocessFunnorm(rgset_ibd, sex=sampleinfo$sex)
ibd_beta<-getBeta(MSet.illumina)





#' ###  beta distribution before and after normalization
set.seed(1)

#extract beta values for plotting
ibd_beta_raw<-getBeta(rgset_ibd)
#remove 6 from raw for plot
ibd_beta_raw<-ibd_beta_raw[,which(colnames(ibd_beta_raw)%in%sampleinfo$array.id)]
identical(colnames(ibd_beta_raw),sampleinfo$array.id)

## for speed with only plot a sample of probes
Beta_melted<- melt(ibd_beta[sample(nrow(ibd_beta),(nrow(ibd_beta)/subset)),])
Beta_melted_raw<- melt(ibd_beta_raw[sample(nrow(ibd_beta),(nrow(ibd_beta)/subset)),])

#remove NAs before plotting (otherwise get many non-inifnite warnings)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
Beta_Plot_raw<-Beta_melted_raw[which(!(is.na(Beta_melted_raw$value))),]

#add meta
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,sampleinfo, by.x="ID", by.y="array.id")
colnames(Beta_Plot_raw)<-c("CpG","ID","Beta")
Beta_Plot_raw<-merge(Beta_Plot_raw,sampleinfo, by.x="ID", by.y="array.id")


ggplot(Beta_Plot, aes(Beta, group=as.character(pdata.id), color=as.character(sample.site)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")

ggsave(here("DNAm/figs","beta_distribution_EPIC.pdf"), w=10, h=5)
ggsave(here("DNAm/figs/jpeg","beta_distribution_EPIC.jpeg"), w=10, h=5)

ggplot(Beta_Plot_raw, aes(Beta, group=as.character(pdata.id), color=as.character(sample.site)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm Beta Value")

ggsave(here("DNAm/figs","beta_distribution_raw_EPIC.pdf"), w=10, h=5)
ggsave(here("DNAm/figs/jpeg","beta_distribution_raw_EPIC.jpeg"), w=10, h=5)





#' ### Clustering By Any Meta Data variable 
#' Similar to the 450K samples clustered by sample site of origin when clustering on the whole array. There is some mixing of TI and SC but that was also seen on the 450K.

# remove rows with NAs
ibdBetas_cluster<-ibd_beta[complete.cases(ibd_beta),]

d <- dist(t(ibdBetas_cluster))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(sampleinfo$sample.site), cex=1.5)

pdf(here("DNAm/figs","cluster_wholeEPIC.pdf"), width=30)
myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(sampleinfo$sample.site), cex=1.5)
dev.off()

hclusters <- cutree(hc, h=78)
table(true=sampleinfo$sample.site, cluster=hclusters)

#' Genotyping Probes
#' Unfortunately when clustering on the 64 SNPs on the EPIC the 4 samples from case 50 do not cluster together. The rescope samples from this individual more closely cluster with samples from individual 206. This can be confirmed later with the genotyping data and check with Fliss to confirm sample IDs.
SNPs <- getSnpBeta(rgset_ibd)
SNPs<-SNPs[complete.cases(SNPs),]# 65 cause one was all NA
SNPs_EPIC<-SNPs[complete.cases(SNPs),]# for saving
save(SNPs_450K, SNPs_EPIC, file=paste(here("DNAm/data"),"/SNPS_on_DNAm.RData",sep=""))# saved to compared to samples ids from omni express exome

SNPs<-SNPs[,which(colnames(SNPs)%in%sampleinfo$array.id)]
identical(colnames(SNPs),sampleinfo$array.id)

d <- dist(t(SNPs))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(as.character(sampleinfo$case.no)), cex=1.5)

pdf(here("DNAm/figs","cluster_snps_EPIC.pdf"), width=30)
  myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(as.character(sampleinfo$case.no)), cex=1.5)
dev.off()

hclusters <- cutree(hc, h=2)
ggplot(melt(table(true=sampleinfo$case.no, cluster=hclusters)), aes(as.character(true), cluster, fill=as.factor(value)))+geom_tile(color='black')+
  scale_fill_manual(values=c('white','#9ecae1','#6baed6','#4292c6','#2171b5','#084594'))
tbl<-table(true=sampleinfo$case.no, cluster=hclusters)
sapply(1:nrow(tbl), function(x) length(which(tbl[x,]!=0)))


#' by sex too
#' When clustering on the sex chromosomes, again sample 50 rescopes, cluster incorrectly. Individual 50 is labelled male but the rescopes cluster with the female samples, 206 is female and therefore this may be more evidence the 50 rescopes could be 206 rescope. Two other samples cluster with the incorrect sex, samples 214 and 65. This pattern in clear on the X chromosome (19,090 CpGs), but less on the Y (537 CpGs).
anno_EPIC<-read.csv(here("../data","MethylationEPIC_v-1-0_B4.csv"), skip=7)
ibd_beta<-ibd_beta[which(rownames(ibd_beta)%in%anno_EPIC$IlmnID),]
anno_EPIC<-anno_EPIC[match(rownames(ibd_beta),anno_EPIC$IlmnID),]

ibd_sex<-ibd_beta[which(anno_EPIC$CHR%in%c('X','Y')),]
d <- dist(t(ibd_sex))
hc <- hclust(d, method = "complete") #single, complete, average, ward

myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(as.character(sampleinfo$case.no)), cex=1.5)

pdf(here("DNAm/figs","cluster_sex_EPIC.pdf"), width=30)
myplclust(hc, labels=sampleinfo$sample_ID, lab.col=as.fumeric(sampleinfo$sex), cex=1.5)
dev.off()

#' #### Another check on sex mixups
sampleinfo$sex_mix<-sapply(1:nrow(sampleinfo), function(x) if(sampleinfo$case.no[x] %in% c("65","50","214")){sampleinfo$case.no[x]}else{sampleinfo$sex[x]})

# beta plot y chr
ibd_Y<-ibd_beta[which(anno_EPIC$CHR%in%c('Y')),]#537
Beta_melted<- melt(ibd_Y)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,sampleinfo, by.x="ID", by.y="array.id")

ggplot(Beta_Plot, aes(Beta, group=as.character(pdata.id), color=as.character(sex_mix)))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")+scale_color_manual(values=c("#a50f15","#ef3b2c","#08519c", "#9ecae1","#fc9272"), name="Sex")
ggsave(here("DNAm/figs","beta_distribution_EPIC_Y.pdf"), w=10, h=5)
ggsave(here("DNAm/figs/jpeg","beta_distribution_EPIC_Y.jpeg"), w=10, h=5)


# beta plot x chr
ibd_X<-ibd_beta[which(anno_EPIC$CHR%in%c('X')),]#19090
Beta_melted<- melt(ibd_X)
Beta_Plot<-Beta_melted[which(!(is.na(Beta_melted$value))),]
colnames(Beta_Plot)<-c("CpG","ID","Beta")
Beta_Plot<-merge(Beta_Plot,sampleinfo, by.x="ID", by.y="array.id")

ggplot(Beta_Plot, aes(Beta, group=as.character(pdata.id), color=as.character(sex_mix)))+
  geom_density()+theme_bw()+xlab("DNAm Beta Value")+scale_color_manual(values=c("#a50f15","#ef3b2c","#08519c", "#9ecae1","#fc9272"), name="Sex")
ggsave(here("DNAm/figs","beta_distribution_EPIC_X.pdf"), w=10, h=5)
ggsave(here("DNAm/figs/jpeg","beta_distribution_EPIC_X.jpeg"), w=10, h=5)

#' ## Sample mixups
#' From Fliss
#' 1. Sample 65: should indeed be labelled as male
#' 2. Sample 214: should be labelled as female (sorry about these 2 – this is the problem with manual record keeping…)
#' 3. Sample 102: should have 11 as the age at diagnosis (not 1) – again, I imagine this will have been a typo at some stage

sampleinfo[grep("65", sampleinfo$case.no),]$sex<-"M"
sampleinfo[grep("214", sampleinfo$case.no),]$sex<-"F"
sampleinfo[grep("102", sampleinfo$case.no),]$age<-"11"

write.table(sampleinfo[,1:16], file=paste(here(path),"/Zerbino_EPICMethylationArray_SampleInfo_RE_mislabelled_fixed.txt",sep=""),sep="\t")




#' ## probe filtering detection pvalue and backgrund
#' SNP probes already filtered by overlapping with EPIC annotation. Sex chromosomes have a very different distribution between sexes so they will be removed prior to normalization. Probes shown to hybridze to multiple parts of the genome will be filtered. Probes with a known SNP at the C of G of the queried CpG will be filtered. Probes with a high number of NA resulting from low beadcount (in 5% of samples) will be removed.  Probes with a high detection p value (>0.05) in 1% of samples will be removed. The number of filtered probes at each stage is shown in the plot below. On the EPIC a very high number of probes were removed for high detection p value. Likely a result of the one array with widespread high detection p values. Removing the 6 bad samples with high detecttion p values on the one array resulted in 30,000 less probes being filtered for high detection p values. 
ibd_beta <- ibd_beta[!grepl("rs",rownames(ibd_beta)), ]
dim(ibd_beta) # probes = 865859, n = 96

#' ## Sex CHR
anno_EPIC<-anno_EPIC[anno_EPIC$IlmnID%in%rownames(ibd_beta),]
ibd_beta <- ibd_beta[!anno_EPIC$CHR%in%c("X", "Y"), ]
dim(ibd_beta) # probes = 846232, n = 96, 19627 filtered
filt_sex<-nrow(ibd_beta)


#' ##Cross hybridizing probes
# 'https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1066-1
#'"43,254 cross-reactive probes with ≥ 47 bp homology with an off-target site, of which 15,782 (36.5 %) are new to the EPIC platform"
#' They include this annotated list in their supplement.
#' wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM2_ESM.csv

cross_reactive<-read.csv(here("../data","13059_2016_1066_MOESM2_ESM.csv"), stringsAsFactors = F)
ibd_beta<-ibd_beta[which(!(rownames(ibd_beta)%in%cross_reactive$PROBE)),]
dim(ibd_beta) # probes = 814600, n = 96, 31632 filtered
filt_cross<-nrow(ibd_beta)

#' ## Polymorphic probes
#'For polymorphic probes I will The Pidsley annotation aswell for "Probes overlapping genetic variants at targeted CpG sites."
#'and "Probes overlapping genetic variants at single base extension sites for Infinium Type I probes"
#'but NOT "Probes with genetic variants overlapping the body of the probe: 48 base pairs for Infinium Type I probes and 49 base pairs for Infinium Type II probes."

#'wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM4_ESM.csv
#'wget https://static-content.springer.com/esm/art%3A10.1186%2Fs13059-016-1066-1/MediaObjects/13059_2016_1066_MOESM5_ESM.csv

polymorphic<-read.csv(here("../data","13059_2016_1066_MOESM4_ESM.csv"),stringsAsFactors = F)
length(unique(polymorphic$PROBE))
baseext<-read.csv(here("../data","13059_2016_1066_MOESM5_ESM.csv"),stringsAsFactors = F)
length(unique(baseext$PROBE))
ibd_beta<-ibd_beta[which(!(rownames(ibd_beta)%in%c(polymorphic$PROBE, baseext$PROBE))),]
dim(ibd_beta) # probes = 802848, n = 96, 11752 filtered
filt_poly<-nrow(ibd_beta)











#' ## probe filtering detection pvalue and backgrund

#' Remove probes with high NA count
na_count_probe <-sapply(1:nrow(ibd_beta), function(y) length(which(is.na(ibd_beta[y,]))))
na_count_probe_good<-which(na_count_probe<(ncol(ibd_beta)*0.05))
ibd_beta<-ibd_beta[na_count_probe_good,]
dim(ibd_beta)# probes = 802848, n = 96, 0 filtered
filt_bead<-nrow(ibd_beta)


#' detection pval
detP <- detectionP(rgset_ibd)
detP<-detP[,which(colnames(detP)%in%sampleinfo$array.id)]
identical(colnames(detP),sampleinfo$array.id)

failed <- detP>0.05
bad_det_p<-names(which(rowMeans(failed)>0.01))
bad_det_psamp<-names(which(colMeans(failed)>0.01))

ibd_beta<-ibd_beta[which(!(rownames(ibd_beta)%in%bad_det_p)),]
ibd_beta<-ibd_beta[,which(!(colnames(ibd_beta)%in%bad_det_psamp))]

dim(ibd_beta)# probes = 767439, n = 96,  35409 filtered
filt_detp<-nrow(ibd_beta)


#' these numbers match pfilter when run on raw rgset data
#' 0 samples having 1 % of sites with a detection p-value greater than 0.05 were removed
#' 0 sites were removed as beadcount <3 in 5 % of samples
#' 35409 sites having 1 % of samples with a detection p-value greater than 0.05 were removed

ibd_beta_epic<-ibd_beta

#' remove 50 rescopes
RE50<-sampleinfo$array.id[which(sampleinfo$case.no=="50" & sampleinfo$passage.or.rescope.no=="RE1")]
ibd_beta_epic<-ibd_beta_epic[,which(!(colnames(ibd_beta_epic)%in%RE50))]
save(ibd_beta_epic, ibd_beta_450K, file=paste(here("../data"),"/ibd_beta_botharrays_funnorm_filtered.RData",sep=""))



#' ### Probe attrition plot
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

ggsave(here("DNAm/figs","probe_attrition_EPIC.pdf"), width = 8, height = 3)
ggsave(here("DNAm/figs/jpeg","probe_attrition_EPIC.jpeg"), width = 8, height = 3)



#' ## PCA
#'To look at overall contributors to EPIC variation PCA was run and PC loadings were associated to meta data variable. Generally the sample site of the epithelial cells is the primary driver of variation, followed by the age and case no. (ie patient).  Interestingly inflammation and sample site are overlapping in PC1.
pca_res <- prcomp(t(ibd_beta))

Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)

#Restructure meta
sampleinfo$sentrix_ID<-as.factor(sampleinfo$sentrix_ID)
sampleinfo$case.no<-as.factor(sampleinfo$case.no)

meta_categorical <- sampleinfo[, c(1,4,5,10,11,18)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(sampleinfo[, c(12)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Case No.", "Diagnosis", "Sample Site","Inflammation","Sex","Sentrix ID")
colnames(meta_continuous) <- c("Age")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7))

ggsave(here("DNAm/figs","heat_scree_EPIC.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7)),width = 9, height = 6)
ggsave(here("DNAm/figs/jpeg","heat_scree_EPIC.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 2.5, 2.7)),width = 9, height = 6)

## PC vs PC plot
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, sampleinfo, by="array.id")

ggplot(Loadings_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=2, color="black")+theme_bw()
ggsave(here("DNAm/figs","PC1_PC2_EPIC.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_EPIC.jpeg"), width = 7.5, height = 6)








#'## R Session Info
sessionInfo()



