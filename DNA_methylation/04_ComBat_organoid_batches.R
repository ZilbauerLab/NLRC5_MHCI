

#'---
#'title: Combat array type but no deconvolution
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---


#' ### Load Libraries
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


#' ### Load Functions
source(here("general_functions","00_pretty_plots.R"))
source(here("general_functions","00_Heat_scree_plot_generic.R"))


#' ## Load first organoid batch
load(here("../data/ibd_beta_organoids.RData"))
epic.organoid$passage.or.rescope.no_numeric<-as.factor(epic.organoid$passage.or.rescope.no)
levels(epic.organoid$passage.or.rescope.no_numeric)<-c(1,10,11,14,16,2,3,4,6,7,8,2,4)
epic.organoid$passage.or.rescope.no_numeric<-as.numeric(as.character(epic.organoid$passage.or.rescope.no_numeric))
epic.organoid$Sample_Name<-paste(epic.organoid$case.no,epic.organoid$sample.site, epic.organoid$passage.or.rescope.no, sep=" ")



#'### load new batch organoids
load(file=here("DNAm/data/newbatch","newbatch_betas_normalized_aftersampleQC.RData"))

#' ## Combine data
head(newbatch_epic.organoid)
newbatch_epic.organoid$array.type<-"EPIC"
newbatch_epic.organoid$sample_ID<-paste(newbatch_epic.organoid$individual,"_",newbatch_epic.organoid$Segment, sep="")
newbatch_epic.organoid$batch<-"new"
newbatch_epic.organoid<-newbatch_epic.organoid[c("individual","array.id","diagnosis_simplified","Segment",
                                                 "Gender","Age","array.type","Sentrix_Position",
                                                 "sample_ID","Sentrix_ID", "det_pval","passage","Sample_Name","Wnt.type",
                                                 "condition","Biobank.Rachel.Replicates","batch")]



head(epic.organoid)
epic.organoid<-epic.organoid[,which(!(colnames(epic.organoid)%in%c("age.group","sample.type","sampling.time.point","passage.or.rescope.no", "fraction","differentiated",
                                                                   "knockout","aza.treatment","inflammation", "ae.4957","ae.5463","id.4957","gpaed.id","tripp.id","other.id","pdata.id",
                                                                   "genotype","notes","long","duplicated","sentrix.id","array.id.path")))]
epic.organoid$Wnt.type<-"WCM"
epic.organoid$condition<-NA
epic.organoid$Biobank.Rachel.Replicates<-"old_Biobank"
epic.organoid$batch<-"old"

colnames(newbatch_epic.organoid)<-colnames(epic.organoid)

newbatch_epic.organoid$sample.site<-as.character(newbatch_epic.organoid$sample.site)
newbatch_epic.organoid$sentrix_ID<-as.character(newbatch_epic.organoid$sentrix_ID)
newbatch_epic.organoid$age<-as.numeric(newbatch_epic.organoid$age)


### Thrid batch organoids
load(file=here("DNAm/data/batch_three_raw","thirdbatch_betas_normalized.RData"))
thirdbatch_epic.organoid$Segment<-as.factor(thirdbatch_epic.organoid$Segment)
levels(thirdbatch_epic.organoid$Segment)<-c("NA","DUO","SC","TI")

## tidier AZA
thirdbatch_epic.organoid$AZA_patient_treated<-thirdbatch_epic.organoid$AZA

## passage numeric
thirdbatch_epic.organoid$passage.or.rescope.no_numeric<-as.factor(thirdbatch_epic.organoid$passage.or.rescope.no)
levels(thirdbatch_epic.organoid$passage.or.rescope.no_numeric)<-c(1,10,12,2,3,4,6,7,8,9)
thirdbatch_epic.organoid$passage.or.rescope.no_numeric<-as.numeric(as.character(thirdbatch_epic.organoid$passage.or.rescope.no_numeric))

## match meta data
thirdbatch_epic.organoid$Wnt.type<-"WCM"
thirdbatch_epic.organoid$batch<-"thrid"

thirdbatch_epic.organoid<-thirdbatch_epic.organoid[c("case.no","array.id","diagnosis","tissue",
                                                     "sex","age","array.type","sentrix.pos",
                                                     "sample.id.3","sentrix.id", "det_pval","passage.or.rescope.no_numeric","Sample_Name","Wnt.type",
                                                     "treatment","notes", "batch")]
colnames(thirdbatch_epic.organoid)<-colnames(epic.organoid)
thirdbatch_epic.organoid$Sample_Name<-gsub("_"," ",thirdbatch_epic.organoid$Sample_Name)



epic.organoid_combined<-rbind(epic.organoid,newbatch_epic.organoid,thirdbatch_epic.organoid)
nrow(epic.organoid_combined)# 455 samples

nrow(ibd_EPIC_organoid_beta)#800383
nrow(newbatch_organoid_beta)#801066
nrow(thirdbatch_organoid_beta)#799126

ibd_EPIC_organoid_beta<-ibd_EPIC_organoid_beta[which(rownames(ibd_EPIC_organoid_beta)%in%rownames(newbatch_organoid_beta)),]
newbatch_organoid_beta<-newbatch_organoid_beta[which(rownames(newbatch_organoid_beta)%in%rownames(ibd_EPIC_organoid_beta)),]
nrow(ibd_EPIC_organoid_beta)#799922
nrow(newbatch_organoid_beta)#799922

identical(rownames(ibd_EPIC_organoid_beta), rownames(newbatch_organoid_beta))

organoid_beta_combined<-cbind(ibd_EPIC_organoid_beta,newbatch_organoid_beta)
dim(organoid_beta_combined)

thirdbatch_organoid_beta<-thirdbatch_organoid_beta[which(rownames(thirdbatch_organoid_beta)%in%rownames(organoid_beta_combined)),]
organoid_beta_combined<-organoid_beta_combined[which(rownames(organoid_beta_combined)%in%rownames(thirdbatch_organoid_beta)),]
nrow(organoid_beta_combined)#797564
nrow(thirdbatch_organoid_beta)#797564

identical(rownames(organoid_beta_combined), rownames(thirdbatch_organoid_beta))

organoid_beta_combined<-cbind(organoid_beta_combined,thirdbatch_organoid_beta)
dim(organoid_beta_combined)#797564    455

identical(colnames(organoid_beta_combined), epic.organoid_combined$array.id)





#' #### replicates cross batches
intersect(newbatch_epic.organoid$Sample_Name, epic.organoid$Sample_Name)
intersect(newbatch_epic.organoid$Sample_Name, thirdbatch_epic.organoid$Sample_Name)
intersect(epic.organoid$Sample_Name, thirdbatch_epic.organoid$Sample_Name)


#' ## PCA for batch effect
pca_res <- prcomp(t(organoid_beta_combined))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)
print(Importance[1:10])

epic.organoid_combined$sentrix_ID<-as.factor(as.character(epic.organoid_combined$sentrix_ID))

meta_categorical <- epic.organoid_combined[, c(1,3,4,5,10,14,17)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(epic.organoid_combined[, c(6, 12)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Individual","Diagnosis","Segment","Gender","Sentrix ID","WNT Type","Batch")
colnames(meta_continuous) <- c( "Age","Passage")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8))

ggsave(here("DNAm/figs/threebatch_heat_scree_organoids_combined.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)
ggsave(here("DNAm/figs/jpeg","threebatch_heat_scree_organoids_combined.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)


## PC vs PC plot
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, epic.organoid_combined, by="array.id")

ggplot(Loadings_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (26%)")+ylab("PC2 (10%)")+th+theme(axis.text = element_text(size=12),
                                              axis.title = element_text(size=14),
                                              plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_sampsite
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC12_site.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC12_site.jpeg"), width = 5, height = 4)


ggplot(Loadings_meta, aes(PC1, PC2, fill=condition))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (26%)")+ylab("PC2 (10%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC12_condition.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC12_condition.jpeg"), width = 5, height = 4)

ggplot(Loadings_meta, aes(PC1, PC2, fill=batch))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (26%)")+ylab("PC2 (10%)")+th+theme(axis.text = element_text(size=12),
                                              axis.title = element_text(size=14),
                                              plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=c("#A5B5BF","#402281","#368556"))
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC12_batch.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC12_batch.jpeg"), width = 5, height = 4)


ggplot(Loadings_meta, aes(PC2, PC3, fill=sample.site))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC2 (10%)")+ylab("PC3 (7%)")+th+theme(axis.text = element_text(size=12),
                                             axis.title = element_text(size=14),
                                             plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_sampsite
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC23_site.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC23_site.jpeg"), width = 5, height = 4)

ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage.or.rescope.no_numeric)))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC2 (10%)")+ylab("PC3 (9%)")+th+theme(axis.text = element_text(size=12),
                                             axis.title = element_text(size=14),
                                             plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=pass_col, name="Passage")
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC23_passage.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC23_passage.jpeg"), width = 5, height = 4)


ggplot(Loadings_meta, aes(PC3, PC4, fill=as.factor(passage.or.rescope.no_numeric)))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC3 (9%)")+ylab("PC4 (7%)")+th+theme(axis.text = element_text(size=12),
                                             axis.title = element_text(size=14),
                                             plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=pass_col, name="Passage")
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC34_passage.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC34_passage.jpeg"), width = 5, height = 4)

ggplot(Loadings_meta, aes(PC3, PC4, fill=batch))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC3 (9%)")+ylab("PC4 (7%)")+th+theme(axis.text = element_text(size=12),
                                             axis.title = element_text(size=14),
                                             plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=c("#A5B5BF","#402281","#368556"))
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC34_batch.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC34_batch.jpeg"), width = 5, height = 4)




###############################
#' ## Combat Batch
###############################
# # impute 0 and 1
organoid_beta_combined[organoid_beta_combined==0]<-0.01
organoid_beta_combined[organoid_beta_combined==1]<-0.99

# impute NA
imputeMedianv3<-function(x) apply(x, 1, function(x){x[is.na(x)]<-median(x, na.rm=T); x}) #impute with row mean
organoid_beta_combined<-t(imputeMedianv3(organoid_beta_combined))



Mval<-function(beta) log2(beta/(1-beta))
edata = apply(organoid_beta_combined, 1, Mval) # need mvalues for combat
edata = as.data.frame(edata)
edata = t(edata)

# how many cross batch replicates
table(epic.organoid_combined$Sample_Name)[which(table(epic.organoid_combined$Sample_Name)>1)]
epic.organoid_combined[which(epic.organoid_combined$Sample_Name%in%names(table(epic.organoid_combined$Sample_Name)[which(table(epic.organoid_combined$Sample_Name)>1)])),c("Sample_Name","Biobank.Rachel.Replicates", "batch")]

mod = model.matrix(~as.factor(Sample_Name), data=epic.organoid_combined)
batch = epic.organoid_combined$batch
combat_organoids_mval = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE)

#Back to betas
betas<-function(M) 2^M/((2^M)+1)
combat_organoid_Beta = apply(combat_organoids_mval, 1, betas) # need mvalues for combat
combat_organoid_Beta = as.data.frame(combat_organoid_Beta)
combat_organoid_Beta = t(combat_organoid_Beta)
combat_organoid_Beta<-as.data.frame(combat_organoid_Beta)

combat_organoid_Beta<-t(imputeMedianv3(combat_organoid_Beta))
save(combat_organoid_Beta,epic.organoid_combined, file=paste(here("DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))







load(file=paste(here("DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))
## fix 202 and 203
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605870061_R03C01")]<-"hold"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605880096_R04C01")]<-"205605870061_R03C01"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="hold")]<-"205605880096_R04C01"
### load update diagnosis
# epic.samples<-read.table(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","AllEPICOrganoid_MHCIannotated_UpdatedSampleInfo_31Oct22.txt"), header=T, sep="\t")
epic.samples<-read.table(here("DNAm/data","AllEPICOrganoid_MHCIannotated_UpdatedSampleInfo_31Oct22.txt"), header=T, sep="\t")
epic.organoid_combined_update<-merge(epic.organoid_combined[,c("array.id" , "sample.site" ,"array.type","sentrix.pos" ,
                                                               "sample_ID","sentrix_ID","det_pval","passage.or.rescope.no_numeric", "Sample_Name",
                                                               "Wnt.type","condition","Biobank.Rachel.Replicates","batch")],
                                     epic.samples[,c("case.no", "array.id","diagnosis", "sex", "age","sample.id","numeric.passage",
                                                     "Disease.Severity","DUO.Inflammation", "TI.Inflammation", "SC.Inflammation",
                                                     "Biologics", "Surgery","AZA", "Perianal.Disease")],
                                     by="array.id")



#' ## PCA for batch effect
pca_res <- prcomp(t(combat_organoid_Beta))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)
print(Importance[1:10])

epic.organoid_combined$sentrix_ID<-as.factor(as.character(epic.organoid_combined$sentrix_ID))

meta_categorical <- epic.organoid_combined[, c(1,3,4,5,10,14,17)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(epic.organoid_combined[, c(6, 12)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Individual","Diagnosis","Segment","Gender","Sentrix ID","WNT Type","Batch")
colnames(meta_continuous) <- c( "Age","Passage")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8))

ggsave(here("DNAm/figs/threebatch_heat_scree_organoids_combined_combat.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)
ggsave(here("DNAm/figs/jpeg","threebatch_heat_scree_organoids_combined_combat.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)


## PC vs PC plot
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, epic.organoid_combined_update, by="array.id")

ggplot(Loadings_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (26%)")+ylab("PC2 (9%)")+th+theme(axis.text = element_text(size=12),
                                              axis.title = element_text(size=14),
                                              plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_sampsite
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC12_site_combat.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC12_site_combat.jpeg"), width = 5, height = 4)

ggplot(Loadings_meta, aes(PC1, PC2, fill=diagnosis))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (26%)")+ylab("PC2 (9%)")+th+theme(axis.text = element_text(size=12),
                                              axis.title = element_text(size=14),
                                              plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_diagnosis
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC12_diagnosis_combat.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC12_diagnosis_combat.jpeg"), width = 5, height = 4)

ggplot(Loadings_meta, aes(PC1, PC2, fill=diagnosis))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (26%)")+ylab("PC2 (9%)")+th+theme(axis.text = element_text(size=12),
                                              axis.title = element_text(size=14),
                                              plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_diagnosis
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC12_diagnosis_combat.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC12_diagnosis_combat.jpeg"), width = 5, height = 4)


ggplot(Loadings_meta, aes(PC2, PC3, fill=sample.site))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC2 (9%)")+ylab("PC3 (9%)")+th+theme(axis.text = element_text(size=12),
                                             axis.title = element_text(size=14),
                                             plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_sampsite
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC23_site_combat.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC23_site_combat.jpeg"), width = 5, height = 4)

ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage.or.rescope.no_numeric)))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC2 (9%)")+ylab("PC3 (9%)")+th+theme(axis.text = element_text(size=12),
                                             axis.title = element_text(size=14),
                                             plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=pass_col, name="Passage")
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC23_passage_combat.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC23_passage_combat.jpeg"), width = 5, height = 4)


ggplot(Loadings_meta, aes(PC3, PC4, fill=as.factor(passage.or.rescope.no_numeric)))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC3 (9%)")+ylab("PC4 (7%)")+th+theme(axis.text = element_text(size=12),
                                             axis.title = element_text(size=14),
                                             plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=pass_col, name="Passage")
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC34_passage_combat.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC34_passage_combat.jpeg"), width = 5, height = 4)

ggplot(Loadings_meta, aes(PC3, PC4, fill=batch))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC3 (9%)")+ylab("PC4 (7%)")+th+theme(axis.text = element_text(size=12),
                                             axis.title = element_text(size=14),
                                             plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=c("#A5B5BF","#402281","#368556"))
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC34_batch_combat.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC34_batch_combat.jpeg"), width = 5, height = 4)





###########
#' ## PCA for batch effect no Duo
###########
# epic.organoid_combined_SCTI<-epic.organoid_combined_update[which(epic.organoid_combined_update$sample.site%in%c("SC","TI")),]
combat_organoid_Beta_SCTI<-combat_organoid_Beta[,which(colnames(combat_organoid_Beta)%in%epic.organoid_combined_SCTI$array.id)]
epic.organoid_combined_SCTI<-epic.organoid_combined_SCTI[match(colnames(combat_organoid_Beta_SCTI),epic.organoid_combined_SCTI$array.id),]
identical(colnames(combat_organoid_Beta_SCTI), epic.organoid_combined_SCTI$array.id)

pca_res <- prcomp(t(combat_organoid_Beta_SCTI))
Loadings<-as.data.frame(pca_res$x)
vars <- pca_res$sdev^2
Importance<-vars/sum(vars)
print(Importance[1:10])


epic.organoid_combined_SCTI$sentrix_ID<-as.factor(as.character(epic.organoid_combined_SCTI$sentrix_ID))

meta_categorical <- epic.organoid_combined_SCTI[, c(1,3,4,5,10,14,17)]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(epic.organoid_combined_SCTI[, c(6, 12)] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Individual","Diagnosis","Segment","Gender","Sentrix ID","WNT Type","Batch")
colnames(meta_continuous) <- c( "Age","Passage")

ord<-1:length(c(colnames(meta_categorical),colnames(meta_continuous)))
# how far do you want the plot to go?
PCs_to_view<-10

suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8))

ggsave(here("DNAm/figs/threebatch_heat_scree_organoids_combined_combat_noduo.pdf"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)
ggsave(here("DNAm/figs/jpeg","threebatch_heat_scree_organoids_combined_combat_noduo.jpeg"), suppressWarnings(heat_scree_plot(Loadings, Importance, 3.3, 1.8)),width = 9, height = 6)


## PC vs PC plot
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, epic.organoid_combined_SCTI, by="array.id")

ggplot(Loadings_meta, aes(PC2, PC3, fill=as.factor(passage.or.rescope.no_numeric)))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC2 (11%)")+ylab("PC3 (4%)")+th+theme(axis.text = element_text(size=12),
                                             axis.title = element_text(size=14),
                                             plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=pass_col, name="Passage")


ggplot(Loadings_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (27%)")+ylab("PC2 (10%)")+th+theme(axis.text = element_text(size=12),
                                             axis.title = element_text(size=14),
                                             plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_sampsite
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC12_samplesite_combat_SCTI.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC12_samplesite_combat_SCTI.jpeg"), width = 5, height = 4)

ggplot(Loadings_meta, aes(PC1, PC2, fill=diagnosis))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC1 (27%)")+ylab("PC2 (10%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  fillscale_diagnosis
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC12_diagnosis_combat_SCTI.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC12_diagnosis_combat_SCTI.jpeg"), width = 5, height = 4)



ggplot(Loadings_meta, aes(PC3, PC4, fill=batch))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC3 (4%)")+ylab("PC4 (3%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=c("#A5B5BF","#402281","#368556"))
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC34_batch_combat_SCTI.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC34_batch_combat_SCTI.jpeg"), width = 5, height = 4)


ggplot(Loadings_meta, aes(PC3, PC4, fill=Wnt.type))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC3 (4%)")+ylab("PC4 (3%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))+
  scale_fill_manual(values=c("darkblue","grey"))
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC34_wnt_combat_SCTI.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC34_wnt_combat_SCTI.jpeg"), width = 5, height = 4)


ggplot(Loadings_meta, aes(PC3, PC4, fill=age))+geom_point(shape=21,size=3, color="black")+theme_bw()+
  xlab("PC3 (4%)")+ylab("PC4 (3%)")+th+theme(axis.text = element_text(size=12),
                                               axis.title = element_text(size=14),
                                               plot.margin = margin(1, 0.1, 1, 1, "cm"))
ggsave(here("DNAm/figs","three_batch_organoids_combined_PC34_batch_age_SCTI.pdf"), width = 5, height = 4)
ggsave(here("DNAm/figs/jpeg","three_batch_organoids_combined_PC34_batch_age_SCTI.jpeg"), width = 5, height = 4)


