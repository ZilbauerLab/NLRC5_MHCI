
#'---
#'title: Estimate the T cell contamination in primary epithelial cells
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---


#' ##  Adapted from the estimateCellCounts function from minfi and ECC7 script by Dr. Meaghan Jones

#' ### Load Libraries
#' require(quadprog) only needed fro more than 2 cell types
library(genefilter)
library(matrixStats)
library(limma)
library(ggplot2)
library(gridExtra)
library(quadprog)
library(sva)
library(pamr)
library(limma)
library(here)

options(stringsAsFactors = FALSE)

source(here("DNAm/scripts/","04_minfi_estcellcounts_functions.R"))
source(here("general_functions","00_pretty_plots.R"))
source(here("general_functions","00_Heat_scree_plot_generic.R"))


                                    
#' ## Load Pediatric IEC data
load(here("../data/ibd_beta_botharrays_funnorm_filtered.RData"))
sampleinfo_EPIC <- read.table(here("../data/raw/EPIC_METHYLATION_DATA/Zerbino_EPICMethylationArray_SampleInfo_RE_mislabelled_fixed.txt"), header=T, sep="\t")
sampleinfo_450K <- read.table(here("../data/raw/K450_METHYLATION_DATA/Zerbino_K450MethylationArray_SampleInfo.txt"), header=T, sep="\t")
sampleinfo_EPIC<-sampleinfo_EPIC[which(sampleinfo_EPIC$array.id%in%colnames(ibd_beta_epic)),]
sampleinfo<-rbind(sampleinfo_450K,sampleinfo_EPIC)
sampleinfo$sample_ID<-sapply(1:nrow(sampleinfo), function(x) paste(sampleinfo$case.no[x],"_",sampleinfo$sample.site[x], sep=""))
sampleinfo$sentrix_ID<-gsub("\\_.*","",sampleinfo$array.id)
ibd_beta_450K<-ibd_beta_450K[which(rownames(ibd_beta_450K)%in%rownames(ibd_beta_epic)),]
ibd_beta_epic<-ibd_beta_epic[which(rownames(ibd_beta_epic)%in%rownames(ibd_beta_450K)),]
ibd_beta_epic<-ibd_beta_epic[match(rownames(ibd_beta_450K),rownames(ibd_beta_epic)),]
identical(rownames(ibd_beta_epic),rownames(ibd_beta_450K))
ibd_combo<-cbind(ibd_beta_450K,ibd_beta_epic)

sampleinfo<-sampleinfo[which(sampleinfo$array.id%in%colnames(ibd_combo)),]
sampleinfo<-sampleinfo[match(colnames(ibd_combo),sampleinfo$array.id),]
identical(sampleinfo$array.id, colnames(ibd_combo))
    
#'#### PC vs PC plot
pca_res <- prcomp(t(ibd_combo))
Loadings<-as.data.frame(pca_res$x)
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, sampleinfo, by="array.id")

ggplot(Loadings_meta, aes(PC1, PC2, fill=sample.site, color=inflammation))+geom_point(shape=21,size=2)+theme_bw()+fillscale_sampsite+scale_color_manual(values=c("black", "white"))
ggsave(here("DNAm/figs","PC1_PC2_botharrays.pdf"), width = 5, height = 3.5)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_botharrays.jpeg"), width = 5, height = 3.5)


#' ## load t cells          
load(here("DNAm/data/public","tcells_combatted_GEO.Rdata")) #197 t cells


#' ## load organoids 
sampleinfo_organoid <- read.table(here("DNAm/data/public/organoid/E-MTAB-4957.sdrf.txt"), header=T, sep="\t")
sampleinfo_organoid<-sampleinfo_organoid[,c(1:15,30:34,38:41)]
load(here("DNAm/data/public/organoid","organoid_beta.Rdata"))
identical(colnames(organoid_beta),sampleinfo_organoid$Assay.Name)
#' dont want fetal organoids
sampleinfo_organoid<-sampleinfo_organoid[which(sampleinfo_organoid$Characteristics.biosource.type.=="organoid" & sampleinfo_organoid$Characteristics.developmental.stage.!="fetal stage"),]
sampleinfo_organoid<-data.frame(sample_ID=sampleinfo_organoid$Assay.Name, CellType="organoid")
organoid_beta<-organoid_beta[,which(colnames(organoid_beta)%in%sampleinfo_organoid$sample_ID)]
identical(colnames(organoid_beta),sampleinfo_organoid$sample_ID)



#' ## Combine organoids and tcells
organoid_beta<-organoid_beta[which(rownames(organoid_beta)%in%rownames(combat_tcell_Beta)),]
combat_tcell_Beta<-combat_tcell_Beta[which(rownames(combat_tcell_Beta)%in%rownames(organoid_beta)),]
organoid_beta<-organoid_beta[match(rownames(combat_tcell_Beta), rownames(organoid_beta)),]
identical(rownames(organoid_beta),rownames(combat_tcell_Beta))
tcell_meta<-data.frame(sample_ID=tcell_sampleinfo$sample_ID, CellType="tcell")

#' ## Define reference data
refMeta<-rbind(sampleinfo_organoid,tcell_meta)
ref<-cbind(organoid_beta, combat_tcell_Beta)

#' ## Define query data (ped IEC)
dataset<-ibd_combo

    
#' ## Define CpG probes in both ref and query    
refprobes <- rownames(ref)  
samp <- dataset[rownames(dataset) %in% refprobes, ] 

#' ## combine ref and query data
commonprobe <- intersect(as.character(rownames(samp)), as.character(refprobes))
comb <- cbind(samp[commonprobe, ], ref[commonprobe, ])
combMeta <- data.frame(sampleNames = c(colnames(samp), colnames(ref)),
                   studyIndex = rep(c("ibd", "reference"), times = c(ncol(samp), ncol(ref))),
                   stringsAsFactors = FALSE)
    
#' ## Quantile Normalizing Data Together
comb.n <- normalizeQuantiles(comb)

    
#' ## Splitting the combined data back to reference and query data after normalization
ref.n <- comb.n[, combMeta$studyIndex == "reference"] 
dim(ref.n)
samp.n <-  comb.n[, combMeta$studyIndex == "ibd"]



#' ## Picking Probes for Composition Estimation
#' ## Estimates from the difference between organoid and blood
compData <- pickCompProbes(ref.n, refMeta, c("organoid","tcell"), 100) # Call the pickCompProbes2 function below to select the probes that can best discern cell types
coefs <- compData$coefEsts
    

#' ## Estimate Composition
#' ## fit a linear model at each CpG
countsUnconstrained <- projectCellType(samp.n[rownames(coefs), ], coefs, nonnegative = TRUE) # Using the weights generated in the last step to estimate the proportion of each cell type
countsConstrained <- projectCellTypeConstrainAlways(samp.n[rownames(coefs), ], coefs, nonnegative = TRUE,lessThanOne = TRUE) # Using the weights generated in the last step to estimate the proportion of each cell type


#' plot
countsUnconstrained<-as.data.frame(countsUnconstrained)
countsUnconstrained$sample_id<-rownames(countsUnconstrained)
sampleinfo_DNAm_countsUnconstrained<-merge(sampleinfo, countsUnconstrained, by.x="array.id",by.y="sample_id")

countsConstrained<-as.data.frame(countsConstrained)
countsConstrained$sample_id<-rownames(countsConstrained)
sampleinfo_DNAm_countsConstrained<-merge(sampleinfo, countsConstrained, by.x="array.id",by.y="sample_id")

save(sampleinfo_DNAm_countsConstrained,sampleinfo_DNAm_countsUnconstrained, file=paste(here("DNAm/data/"),"DNAm_organoid_counts.RData", sep=""))





#load(file=here("DNAm/data","DNAm_organoid_counts.RData"))

#' # check constrained an unconstrained are related
cor.test(sampleinfo_DNAm_countsConstrained$organoid, sampleinfo_DNAm_countsUnconstrained$organoid) # 0.97
cor.test(sampleinfo_DNAm_countsConstrained$tcell, sampleinfo_DNAm_countsUnconstrained$tcell) # 0.98

#' ## Explored difference in composistion with meta data
sampleinfo_DNAm<-sampleinfo_DNAm_countsConstrained
sampleinfo_DNAm$diagnosis_grouped<-as.factor(sampleinfo_DNAm$diagnosis)
levels(sampleinfo_DNAm$diagnosis_grouped)<-c("CD","Control","CD","UC","UC")
sampleinfo_DNAm$diagnosis_grouped<-factor(sampleinfo_DNAm$diagnosis_grouped, levels=c("Control", "CD", "UC"))
sampleinfo_DNAm$sample.site_grouped<-as.factor(sampleinfo_DNAm$sample.site)
levels(sampleinfo_DNAm$sample.site_grouped)<-c("colon","colon","ileum")


tcell<-ggplot(sampleinfo_DNAm, aes(diagnosis_grouped, tcell))+geom_boxplot(fill="white", outlier.shape=NA)+
  geom_point(aes(fill=diagnosis_grouped),size=3,shape=21,position=position_jitter(w=0.2))+
  facet_wrap(~sample.site_grouped)+th+fillscale_diagnosis+theme_bw()+ylim(0,1)+xlab("Diagnosis")+ylab("T Cell Proportion Estimate")+theme(legend.position="none")+th_present

summary(aov(sampleinfo_DNAm$tcell~sampleinfo_DNAm$diagnosis_grouped+sampleinfo_DNAm$sample.site_grouped))

organoid<-ggplot(sampleinfo_DNAm, aes(diagnosis_grouped, organoid))+geom_boxplot(fill="white", outlier.shape=NA)+
  geom_point(aes(fill=diagnosis_grouped),size=3,shape=21,position=position_jitter(w=0.2))+
  facet_wrap(~sample.site_grouped)+th+fillscale_diagnosis+theme_bw()+ylim(0,1)+xlab("Diagnosis")+ylab("Epithelial Proportion Estimate")+theme(legend.position="none")+th_present

summary(aov(sampleinfo_DNAm$organoid~sampleinfo_DNAm$diagnosis_grouped+sampleinfo_DNAm$sample.site_grouped))

grid.arrange(tcell, organoid)

ggsave(here("DNAm/figs","deconvolution_counts.pdf"),grid.arrange(tcell, organoid), width = 7, height = 9)
ggsave(here("DNAm/figs/jpeg","deconvolution_counts.jpeg"),grid.arrange(tcell, organoid), width = 7, height = 9)

sampleinfo_DNAm$sample.site_grouped<-as.factor(sampleinfo_DNAm$sample.site_grouped)
levels(sampleinfo_DNAm$sample.site_grouped)<-c("Sigmoid Colon", "Terminal Ileum")
tcell<-ggplot(sampleinfo_DNAm, aes(diagnosis_grouped, tcell))+geom_boxplot(fill="white", outlier.shape=NA)+
  geom_point(aes(fill=diagnosis_grouped),size=3,shape=21,position=position_jitter(w=0.2))+
  facet_wrap(~sample.site_grouped)+th+fillscale_diagnosis+theme_bw()+ylim(0,1)+xlab("Diagnosis")+ylab("Intraepithelial Lymphocyte\nProportion Estimate")+theme(legend.position="none")+th_present

ggsave(here("DNAm/figs","deconvolution_counts_epithelial.pdf"),tcell, width = 7, height = 4.5)
ggsave(here("DNAm/figs/jpeg","deconvolution_counts_epithelial.jpeg"),tcell, width = 7, height = 4.5)



#' Significance of difference (t test)
tissue<-c("colon", "ileum")
diagnosis<-c("Control", "UC","CD")
combos<-combn(diagnosis, 2)

sig_diff<-do.call(rbind,lapply(1:2, function(x){do.call(rbind,lapply(1:3, function(y){
compare<-sampleinfo_DNAm[which(sampleinfo_DNAm$diagnosis_grouped%in%combos[,y] & sampleinfo_DNAm$sample.site_grouped==tissue[x]),]
data.frame(tissue=tissue[x], diagnosis=paste(combos[,y], collapse=":"),pval=round(t.test(compare$tcell~compare$diagnosis_grouped)$p.value,4))}))}))
sig_diff$padjust<-p.adjust(sig_diff$pval, n=6)
sig_diff



#' ## Composistion consistent in rescopes?
sampleinfo_rescope<-sampleinfo_DNAm[which(sampleinfo_DNAm$passage.or.rescope.no=="RE1"),]
rescope_samples<-unique(sampleinfo_rescope$sample_ID)
sampleinfo_rescope<-sampleinfo_DNAm[which(sampleinfo_DNAm$sample_ID%in%rescope_samples),]

ggplot(sampleinfo_rescope, aes(passage.or.rescope.no, organoid, group=case.no))+
  geom_point(aes(fill=diagnosis_grouped),size=2,shape=21)+geom_line(color="grey")+
  facet_wrap(~sample.site_grouped)+th+fillscale_diagnosis+theme_bw()+ylim(0,1)

ggsave(here("DNAm/figs","deconvolution_rescope_counts.pdf"), width = 7, height = 4)
ggsave(here("DNAm/figs/jpeg","deconvolution_rescope_counts.jpeg"), width = 7, height = 4)





#' ## consistent with histological measure of inflammation
ggplot(sampleinfo_DNAm, aes(inflammation, tcell))+geom_boxplot(fill="white", outlier.shape=NA)+
  geom_point(aes(fill=inflammation),size=2,shape=21,position=position_jitter(w=0.2))+
  facet_wrap(~sample.site_grouped)+th+theme_bw()+ylim(0,1)+xlab("Inflammation")+ylab("T Cell Estimate")+
  scale_fill_manual(values=c("floralwhite","steelblue4"))+th_present+theme(legend.position="none")

t.test(sampleinfo_DNAm$organoid~sampleinfo_DNAm$inflammation)$p.value
t.test(sampleinfo_DNAm[which(sampleinfo_DNAm$sample.site=="SC"),]$organoid~sampleinfo_DNAm[which(sampleinfo_DNAm$sample.site=="SC"),]$inflammation)$p.value
t.test(sampleinfo_DNAm[which(sampleinfo_DNAm$sample.site=="TI"),]$organoid~sampleinfo_DNAm[which(sampleinfo_DNAm$sample.site=="TI"),]$inflammation)$p.value

ggsave(here("DNAm/figs","deconvolution_inf_counts.pdf"), width = 7, height = 4)
ggsave(here("DNAm/figs/jpeg","deconvolution_inf_counts.jpeg"), width = 7, height = 4)


ggplot(sampleinfo_DNAm, aes(inflammation, tcell))+geom_boxplot(fill="white", outlier.shape=NA)+
  geom_point(aes(fill=inflammation),size=3,shape=21,position=position_jitter(w=0.2))+
  facet_wrap(~sample.site_grouped)+th+theme_bw()+ylim(0,1)+xlab("Inflammation")+ylab("Intraepithelial Lymphocyte\nProportion Estimate")+
  scale_fill_manual(values=c("floralwhite","steelblue4"))+th_present+theme(legend.position="none")

ggsave(here("DNAm/figs","deconvolution_counts_epithelial_inflammation.pdf"), width = 7, height = 4.5)
ggsave(here("DNAm/figs/jpeg","deconvolution_counts_epithelial_inflammation.jpeg"), width = 7, height = 4.5)





#' ## cell comp and inf related (in all yes but what about only in IBD?)
sampleinfo_DNAm_IBD<-sampleinfo_DNAm[which(sampleinfo_DNAm$diagnosis!="Control"),]
ggplot(sampleinfo_DNAm_IBD, aes(inflammation, tcell))+geom_boxplot(fill="white", outlier.shape=NA)+
  geom_point(aes(fill=inflammation),size=2,shape=21,position=position_jitter(w=0.2))+
  facet_wrap(~sample.site_grouped)+th+theme_bw()+ylim(0,1)+xlab("Diagnosis")+ylab("T Cell Estimate")

t.test(sampleinfo_DNAm_IBD$organoid~sampleinfo_DNAm_IBD$inflammation)$p.value
t.test(sampleinfo_DNAm_IBD[which(sampleinfo_DNAm_IBD$sample.site=="SC"),]$organoid~sampleinfo_DNAm_IBD[which(sampleinfo_DNAm_IBD$sample.site=="SC"),]$inflammation)$p.value
t.test(sampleinfo_DNAm_IBD[which(sampleinfo_DNAm_IBD$sample.site=="TI"),]$organoid~sampleinfo_DNAm_IBD[which(sampleinfo_DNAm_IBD$sample.site=="TI"),]$inflammation)$p.value

ggsave(here("DNAm/figs","deconvolution_inf_counts_IBD_only.pdf"), width = 7, height = 4)
ggsave(here("DNAm/figs/jpeg","deconvolution_inf_counts_IBD_only.jpeg"), width = 7, height = 4)



#' ## check that composistion associates to PC
pca_res_all <- prcomp(t(ibd_combo))
vars <- pca_res_all$sdev^2
Importance<-vars/sum(vars)
round(Importance, digits=2)

Loadings_all<-as.data.frame(pca_res_all$x)
Loadings_all$array.id<-rownames(Loadings_all)
Loadings_all_meta<-merge(Loadings_all, sampleinfo_DNAm, by="array.id")

ggplot(Loadings_all_meta, aes(PC1, PC2, fill=organoid))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+
  scale_fill_gradientn(colours = c("#0b4899","cornflowerblue","white"), values=c(0,0.5,1),
                       name="Organoid Estimate")+xlab("PC1 (32%)")+ylab("PC2 (20%)")

ggsave(here("DNAm/figs","deconvolution_PCA_counts.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","deconvolution_PCA_counts.jpeg"), width = 7.5, height = 6)



#' # Cell composition adjustment
sampleinfo_DNAm_countsConstrained<-sampleinfo_DNAm_countsConstrained[match(colnames(ibd_combo),sampleinfo_DNAm_countsConstrained$array.id),]
identical(colnames(ibd_combo), sampleinfo_DNAm_countsConstrained$array.id)

#' Impute NAs?
# imputeMedianv3<-function(x) apply(x, 1, function(x){x[is.na(x)]<-median(x, na.rm=T); x}) #impute with row mean
# Blood_beta_imputed<-t(imputeMedianv3(Blood_beta))

avebeta.lm<-apply(ibd_combo, 1, function(x){
  lm(x~sampleinfo_DNAm_countsConstrained$organoid)
})
residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
colnames(residuals)<-colnames(ibd_combo)
ibd_combo_adjusted<-residuals+matrix(apply(ibd_combo, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
save(ibd_combo_adjusted, file=paste(here("DNAm/data/"),"ibd_adjusted_betas.RData", sep=""))




#load("../../data/ibd_adjusted_betas.RData")

#' ## check that composistion does not associates to PC after correction
pca_res_all <- prcomp(t(ibd_combo_adjusted))
vars <- pca_res_all$sdev^2
Importance<-vars/sum(vars)
round(Importance, digits=2)

Loadings_all<-as.data.frame(pca_res_all$x)
Loadings_all_meta<-Loadings_all
Loadings_all_meta$array.id<-rownames(Loadings_all)
Loadings_all_meta<-merge(Loadings_all_meta, sampleinfo_DNAm, by="array.id")

count<-ggplot(Loadings_all_meta, aes(PC1, PC2, fill=organoid))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+
  scale_fill_gradientn(colours = c("#0b4899","cornflowerblue","white"), values=c(0,0.5,1),
                       name="Organoid\nEstimate")+xlab("PC1 (34%)")+ylab("PC2 (11%)")
site<-ggplot(Loadings_all_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+fillscale_sampsite+xlab("PC1 (34%)")+ylab("PC2 (11%)")
diagnosis<-ggplot(Loadings_all_meta, aes(PC1, PC2, fill=diagnosis_grouped))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+fillscale_diagnosis+xlab("PC1 (34%)")+ylab("PC2 (11%)")
batch<-ggplot(Loadings_all_meta, aes(PC1, PC2, fill=array.type))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+xlab("PC1 (34%)")+ylab("PC2 (11%)")

grid.arrange(count,site, diagnosis, batch, ncol=1)

ggsave(here("DNAm/figs","deconvolution_PCA_adjusted.pdf"), grid.arrange(count,site, diagnosis, batch, ncol=1),width = 7.5, height = 24)
ggsave(here("DNAm/figs/jpeg","deconvolution_PCA_adjusted.jpeg"),grid.arrange(count,site, diagnosis, batch, ncol=1), width = 7.5, height = 24)

ggsave(here("DNAm/figs","deconvolution_PCA_counts_adjusted.pdf"), count, width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","deconvolution_PCA_counts_adjusted.jpeg"), count, width = 7.5, height = 6)

ggsave(here("DNAm/figs","deconvolution_PCA_site_adjusted.pdf"), site, width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","deconvolution_PCA_site_adjusted.jpeg"), site, width = 7.5, height = 6)


PC23<-ggplot(Loadings_all_meta, aes(PC2, PC3, fill=array.type))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+xlab("PC2 (11%)")+ylab("PC3 (6%)")
PC34<-ggplot(Loadings_all_meta, aes(PC3, PC4, fill=array.type))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+xlab("PC3 (6%)")+ylab("PC4 (4%)")
ggsave(here("DNAm/figs","PC234_DNAm_adjusted.pdf"), grid.arrange(PC23, PC34), width = 7.5, height = 12)
ggsave(here("DNAm/figs/jpeg","PC234_DNAm_adjusted.jpeg"), grid.arrange(PC23, PC34), width = 7.5, height = 12)


#'## Heat scree plot
#Restructure meta
Loadings_all<-Loadings_all[match(sampleinfo_DNAm$array.id, rownames(Loadings_all)),]

sampleinfo_DNAm$sentrix_ID<-as.factor(sampleinfo_DNAm$sentrix_ID)
sampleinfo_DNAm$case.no<-as.factor(sampleinfo_DNAm$case.no)

meta_categorical <- sampleinfo_DNAm[, c(2,10,11,13,18, 21, 22)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(sampleinfo_DNAm[, c(12,19)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Case No.","Inflammation","Sex","Array Type","Sentrix ID","Diagnosis","Sample Site")
colnames(meta_continuous) <- c("Age","Organoid\nEstimate")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings_all, Importance, 2.5, 2.7))

ggsave(here("DNAm/figs","heat_scree_DNAm_adjusted.pdf"), suppressWarnings(heat_scree_plot(Loadings_all, Importance, 3.35, 1.5)),width = 9, height = 6)
ggsave(here("DNAm/figs/jpeg","heat_DNAm_adjusted.jpeg"), suppressWarnings(heat_scree_plot(Loadings_all, Importance,  3.35, 1.5)),width = 9, height = 6)






#'## Combat array type from adjusted betas
#load("../../data/ibd_adjusted_betas.RData")
sampleinfo_DNAm<-sampleinfo_DNAm[match(colnames(ibd_combo_adjusted),sampleinfo_DNAm$array.id),]
identical(colnames(ibd_combo_adjusted), as.character(sampleinfo_DNAm$array.id))

#' impute 0 and 1
ibd_combo_adjusted[ibd_combo_adjusted==0]<-0.01
ibd_combo_adjusted[ibd_combo_adjusted==1]<-0.99

#' impute NA
imputeMedianv3<-function(x) apply(x, 1, function(x){x[is.na(x)]<-median(x, na.rm=T); x}) #impute with row mean
ibd_combo_adjusted<-t(imputeMedianv3(ibd_combo_adjusted))

Mval<-function(beta) log2(beta/(1-beta))
edata = suppressWarnings(apply(ibd_combo_adjusted, 1, Mval)) # need mvalues for combat)
edata = as.data.frame(edata)
edata = t(edata)



#mod = model.matrix(~as.factor(Tissue), data=tcell_sampleinfo) # can not protect cause confounded
batch = sampleinfo_DNAm$array.type
combat_ibd_mval = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE)

#Back to betas
betas<-function(M) 2^M/((2^M)+1)
combat_ibd_Beta = apply(combat_ibd_mval, 1, betas) # need mvalues for combat
combat_ibd_Beta = as.data.frame(combat_ibd_Beta)
combat_ibd_Beta = t(combat_ibd_Beta)
combat_ibd_Beta<-as.data.frame(combat_ibd_Beta)

combat_ibd_Beta<-t(imputeMedianv3(combat_ibd_Beta))
save(combat_ibd_Beta,sampleinfo_DNAm, file=paste(here("DNAm/data/"),"ibd_adjusted_combatted.Rdata", sep=""))



load(here("DNAm/data","ibd_adjusted_combatted.Rdata"))

#'## Heat scree plot after combat for array type
pca_res_all <- prcomp(t(combat_ibd_Beta))
vars <- pca_res_all$sdev^2
Importance<-vars/sum(vars)
Loadings_all<-as.data.frame(pca_res_all$x)
Loadings_all<-Loadings_all[match(sampleinfo_DNAm$array.id, rownames(Loadings_all)),]

sampleinfo_DNAm$sentrix_ID<-as.factor(sampleinfo_DNAm$sentrix_ID)
sampleinfo_DNAm$case.no<-as.factor(sampleinfo_DNAm$case.no)

meta_categorical <- sampleinfo_DNAm[, c(2,10,11,13,18, 21, 22)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(sampleinfo_DNAm[, c(12,19)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Case No.","Inflammation","Sex","Array Type","Sentrix ID","Diagnosis","Sample Site")
colnames(meta_continuous) <- c("Age","Organoid\nEstimate")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings_all, Importance, 2.5, 2.7))

ggsave(here("DNAm/figs","heat_scree_DNAm_adjusted_combatted.pdf"), suppressWarnings(heat_scree_plot(Loadings_all, Importance, 3.35, 1.5)),width = 9, height = 6)
ggsave(here("DNAm/figs/jpeg","heat_DNAm_adjusted_combatted.jpeg"), suppressWarnings(heat_scree_plot(Loadings_all, Importance,  3.35, 1.5)),width = 9, height = 6)


#' ## array effect gone?
Loadings_all_meta<-Loadings_all
Loadings_all_meta$array.id<-rownames(Loadings_all)
Loadings_all_meta<-merge(Loadings_all_meta, sampleinfo_DNAm, by="array.id")

PC12<-ggplot(Loadings_all_meta, aes(PC1, PC2, fill=array.type))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+xlab("PC1 (38%)")+ylab("PC2 (7%)")
PC23<-ggplot(Loadings_all_meta, aes(PC2, PC3, fill=array.type))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+xlab("PC2 (7%)")+ylab("PC3 (4%)")
PC34<-ggplot(Loadings_all_meta, aes(PC3, PC4, fill=array.type))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+xlab("PC3 (4%)")+ylab("PC4 (43%)")

ggsave(here("DNAm/figs","PC234_DNAm_adjusted.pdf"), grid.arrange(PC12,PC23, PC34), width = 7.5, height = 18)
ggsave(here("DNAm/figs/jpeg","PC234_DNAm_adjusted.jpeg"), grid.arrange(PC12,PC23, PC34), width = 7.5, height = 18)


ggplot(Loadings_all_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=3,color="black")+theme_bw()+th+xlab("PC1 (38%)")+ylab("PC2 (7%)")+fillscale_sampsite

ggsave(here("DNAm/figs","PC1_PC2_combatted.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_combatted.jpeg"), width = 7.5, height = 6)





#'## R Session Info
sessionInfo()