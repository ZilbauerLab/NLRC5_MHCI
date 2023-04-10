#'---
#'title: Explore Variation in Primary samples
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---


#' ### Load Libraries
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(limma)
library(here)
library(gridExtra)


options(stringsAsFactors = FALSE)
options(scipen = 999)


#' ### Load Functions
source(here("general_functions","00_pretty_plots.R"))
source(here("general_functions","00_Heat_scree_plot_generic.R"))


#' ### Load Data
load(here("../data","ibd_beta_botharrays_funnorm_filtered.RData"))




#' #### Combine meta data from both array types 
sampleinfo_EPIC <- read.table(here("../data/raw/EPIC_METHYLATION_DATA","Zerbino_EPICMethylationArray_SampleInfo_RE_mislabelled_fixed.txt"), header=T, sep="\t")
sampleinfo_450K <- read.table(here("../data/raw/K450_METHYLATION_DATA","Zerbino_K450MethylationArray_SampleInfo.txt"), header=T, sep="\t")

sampleinfo_EPIC<-sampleinfo_EPIC[which(sampleinfo_EPIC$array.id%in%colnames(ibd_beta_epic)),]

sampleinfo<-rbind(sampleinfo_450K,sampleinfo_EPIC)
sampleinfo$sample_ID<-sapply(1:nrow(sampleinfo), function(x) paste(sampleinfo$case.no[x],"_",sampleinfo$sample.site[x], sep=""))
sampleinfo$sentrix_ID<-gsub("\\_.*","",sampleinfo$array.id)


#' ### overlapping probes
ibd_beta_450K<-ibd_beta_450K[which(rownames(ibd_beta_450K)%in%rownames(ibd_beta_epic)),]
ibd_beta_epic<-ibd_beta_epic[which(rownames(ibd_beta_epic)%in%rownames(ibd_beta_450K)),]
ibd_beta_epic<-ibd_beta_epic[match(rownames(ibd_beta_450K),rownames(ibd_beta_epic)),]
identical(rownames(ibd_beta_epic),rownames(ibd_beta_450K))
ibd_combo<-cbind(ibd_beta_450K,ibd_beta_epic)

sampleinfo<-sampleinfo[which(sampleinfo$array.id%in%colnames(ibd_combo)),]
sampleinfo<-sampleinfo[match(colnames(ibd_combo),sampleinfo$array.id),]
identical(sampleinfo$array.id, colnames(ibd_combo))




#' ### IBD vs controls in each tissue
sampleinfo$diagnosis_simple<-as.factor(sampleinfo$diagnosis)
levels(sampleinfo$diagnosis_simple)<-c("IBD","Control","IBD","IBD","IBD")

#' More specific grouping
sampleinfo$diagnosis_grouped<-as.factor(sampleinfo$diagnosis)
levels(sampleinfo$diagnosis_grouped)<-c("CD","Control","CD","UC","UC")
sampleinfo$diagnosis_grouped<-factor(sampleinfo$diagnosis_grouped, levels=c("Control", "CD", "UC"))

# More simple grouping of samplesite
sampleinfo$sample.site_grouped<-as.factor(sampleinfo$sample.site)
levels(sampleinfo$sample.site_grouped)<-c("colon","colon","ileum")


#' ### controls only
sampleinfo_controls<-sampleinfo[which(sampleinfo$diagnosis=="Control"),]
ibd_combo_controls<-ibd_combo[,which(colnames(ibd_combo)%in%sampleinfo_controls$array.id)]
identical(colnames(ibd_combo_controls),sampleinfo_controls$array.id)

sampleinfo_ibd<-sampleinfo[which(sampleinfo$diagnosis!="Control"),]
ibd_combo_ibd<-ibd_combo[,which(colnames(ibd_combo)%in%sampleinfo_ibd$array.id)]
identical(colnames(ibd_combo_ibd),sampleinfo_ibd$array.id)






#' ### PCA ibd ctrl all
#' The PCA on both arrays have PC1 and PC2 not orthogonal. I looked into why that might be. First I separated the samples by controls and cases.
pca_res_ctrl <- prcomp(t(ibd_combo_controls))
Loadings_ctrl<-as.data.frame(pca_res_ctrl$x)

pca_res_ibd <- prcomp(t(ibd_combo_ibd))
Loadings_ibd<-as.data.frame(pca_res_ibd$x)

pca_res_all <- prcomp(t(ibd_combo))
Loadings_all<-as.data.frame(pca_res_all$x)


#' PC vs PC plot
Loadings_ctrl$array.id<-rownames(Loadings_ctrl)
Loadings_ctrl_meta<-merge(Loadings_ctrl, sampleinfo_controls, by="array.id")
ctrl<-ggplot(Loadings_ctrl_meta, aes(PC1, PC2, fill=diagnosis_grouped, color=sample.site_grouped))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis+scale_color_manual(values=c("black","white"))

Loadings_ibd$array.id<-rownames(Loadings_ibd)
Loadings_ibd_meta<-merge(Loadings_ibd, sampleinfo_ibd, by="array.id")
ibd<-ggplot(Loadings_ibd_meta, aes(PC1, PC2, fill=diagnosis_grouped, color=sample.site_grouped))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis+scale_color_manual(values=c("black","white"))

Loadings_all$array.id<-rownames(Loadings_all)
Loadings_all_meta<-merge(Loadings_all, sampleinfo, by="array.id")
all<-ggplot(Loadings_all_meta, aes(PC1, PC2, fill=diagnosis_grouped, color=sample.site_grouped))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis+scale_color_manual(values=c("black","white"))

grid.arrange(all,ctrl, ibd)

PC_ctrl_case<-grid.arrange(all,ctrl, ibd)

ggsave(here("DNAm/figs","PC1_PC2_botharrays_split_ctrl_case.pdf"), PC_ctrl_case, width = 7.5, height = 18)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_botharrays_split_ctrl_case.jpeg"), PC_ctrl_case, width = 7.5, height = 18)



#' all diagnosis and sample site for comparison
vars <- pca_res_all$sdev^2
Importance<-vars/sum(vars)

ggplot(Loadings_all_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=3,color="black")+theme_bw()+fillscale_sampsite+th+xlab("PC1 (32%)")+ylab("PC2 (20%)")
ggsave(here("DNAm/figs","PC1_PC2_botharrays_sample_site.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_botharrays_sample_site.jpeg"), width = 7.5, height = 6)

ggplot(Loadings_all_meta, aes(PC1, PC2, fill=diagnosis_grouped, color=sample.site))+geom_point(shape=21,size=3)+theme_bw()+
  fillscale_diagnosis+th+scale_color_manual(values=c("grey40","black","white"), name="Sample Site")

ggsave(here("DNAm/figs","PC1_PC2_botharrays_diagnosis.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_botharrays_diagnosis.jpeg"), width = 7.5, height = 6)


#' sample site simplified and diagnosis
summary(pca_res_all)

levels(Loadings_all_meta$sample.site_grouped)<-c("Colon", "Small\nIntestine")
ggplot(Loadings_all_meta, aes(PC1, PC2, fill=diagnosis_grouped, color=sample.site_grouped))+geom_point(shape=21,size=3)+theme_bw()+
  fillscale_diagnosis+th+scale_color_manual(values=c("black","grey"), name="Sample Site")+xlab("PC1 (32%)")+ylab("PC2 (20%)")

ggsave(here("DNAm/figs","PC1_PC2_botharrays_diagnosis_simpleSampleSite.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_botharrays_diagnosis_simpleSampleSite.jpeg"), width = 7.5, height = 6)



#' ### Density plot wth inflammation (histological measured)
#' So only in the cases are the PCs orthogonal. It is something specific to IBD. However it isn't just histological and endoscope inflammation.
densityinf<-ggplot(Loadings_all_meta, aes(inflammation,PC1, fill=inflammation))+
  geom_violin(alpha=0.5)+geom_boxplot(width=0.1)+th+theme_bw()+
  scale_fill_manual(values=c("lightgrey","#e34a33"))+coord_flip()+theme(plot.margin=unit(c(0.5,1.6,0.5,0.45),"cm"))
scatter<-ggplot(Loadings_all_meta, aes(PC1, PC2, fill=diagnosis_grouped, color=inflammation))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis+scale_color_manual(values=c("black", "white"))

grid.arrange(densityinf ,scatter, heights=c(1.4, 4))
inf_sactter<-grid.arrange(densityinf, scatter, heights=c(1.4, 4))

ggsave(here("DNAm/figs","PC1_PC2_botharrays_cases_inflammation.pdf"),inf_sactter, width = 7.5, height = 8)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_botharrays_cases_inflammation.jpeg"),inf_sactter, width = 7.5, height = 8)



#' ## individual disease
#'  Is the pattern present when looking at UC and CD seperately? 
sampleinfo_UC<-sampleinfo[which(sampleinfo$diagnosis=="UC"),]
ibd_combo_UC<-ibd_combo[,which(colnames(ibd_combo)%in%sampleinfo_UC$array.id)]
identical(colnames(ibd_combo_UC),sampleinfo_UC$array.id)

sampleinfo_CD<-sampleinfo[which(sampleinfo$diagnosis=="CD"),]
ibd_combo_CD<-ibd_combo[,which(colnames(ibd_combo)%in%sampleinfo_CD$array.id)]
identical(colnames(ibd_combo_CD),sampleinfo_CD$array.id)

pca_res_UC <- prcomp(t(ibd_combo_UC))
Loadings_UC<-as.data.frame(pca_res_UC$x)
Loadings_UC$array.id<-rownames(Loadings_UC)
Loadings_UC_meta<-merge(Loadings_UC, sampleinfo_UC, by="array.id")
UC<-ggplot(Loadings_UC_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=2)+theme_bw()+fillscale_sampsite

pca_res_CD <- prcomp(t(ibd_combo_CD))
Loadings_CD<-as.data.frame(pca_res_CD$x)
Loadings_CD$array.id<-rownames(Loadings_CD)
Loadings_CD_meta<-merge(Loadings_CD, sampleinfo_CD, by="array.id")
CD<-ggplot(Loadings_CD_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=2)+theme_bw()+fillscale_sampsite
grid.arrange(CD, UC)



#' ## color by all meta data variables possibly explain variance 
pca_res <- prcomp(t(ibd_combo))
Loadings<-as.data.frame(pca_res$x)
Loadings$array.id<-rownames(Loadings)
Loadings_meta<-merge(Loadings, sampleinfo, by="array.id")

ggplot(Loadings_meta, aes(PC1, PC2, fill=sample.site, color=inflammation))+geom_point(shape=21,size=2)+theme_bw()+fillscale_sampsite+scale_color_manual(values=c("black", "white"))
ggplot(Loadings_meta, aes(PC1, PC2, fill=diagnosis_grouped, color=inflammation))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis+scale_color_manual(values=c("black", "white"))
ggplot(Loadings_meta, aes(PC1, PC2, fill=array.type, color=inflammation))+geom_point(shape=21,size=2)+theme_bw()+scale_color_manual(values=c("black", "white"))
ggplot(Loadings_meta, aes(PC1, PC2, fill=age, color=inflammation))+geom_point(shape=21,size=2)+theme_bw()+scale_color_manual(values=c("black", "white"))
ggplot(Loadings_meta, aes(PC1, PC2, fill=sex, color=inflammation))+geom_point(shape=21,size=2)+theme_bw()+scale_color_manual(values=c("black", "white"))


ggplot(Loadings_meta, aes(PC2, PC3, fill=sample.site, color=inflammation))+geom_point(shape=21,size=2)+theme_bw()+fillscale_sampsite+scale_color_manual(values=c("black", "white"))
ggplot(Loadings_meta, aes(PC2, PC3, fill=diagnosis_grouped, color=inflammation))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis+scale_color_manual(values=c("black", "white"))
ggplot(Loadings_meta, aes(PC2, PC3, fill=array.type, color=diagnosis_simple))+geom_point(shape=21,size=2)+theme_bw()+scale_color_manual(values=c("black", "white"))

ggplot(Loadings_meta, aes(PC1, PC2, fill=as.character(case.no), color=inflammation))+geom_point(shape=21,size=2)+theme_bw()+scale_color_manual(values=c("black", "white"))








#' ### PCA split by array
sampleinfo_450k<-sampleinfo[which(sampleinfo$array.type=="K450"),]
ibd_combo_450k<-ibd_combo[,which(colnames(ibd_combo)%in%sampleinfo_450k$array.id)]
identical(colnames(ibd_combo_450k),sampleinfo_450k$array.id)

sampleinfo_EPIC<-sampleinfo[which(sampleinfo$array.type=="EPIC"),]
ibd_combo_EPIC<-ibd_combo[,which(colnames(ibd_combo)%in%sampleinfo_EPIC$array.id)]
identical(colnames(ibd_combo_EPIC),sampleinfo_EPIC$array.id)


pca_res_450k <- prcomp(t(ibd_combo_450k))
Loadings_450k<-as.data.frame(pca_res_450k$x)

pca_res_EPIC <- prcomp(t(ibd_combo_EPIC))
Loadings_EPIC<-as.data.frame(pca_res_EPIC$x)

## PC vs PC plot
Loadings_450k$array.id<-rownames(Loadings_450k)
Loadings_450k_meta<-merge(Loadings_450k, sampleinfo_450k, by="array.id")
k450<-ggplot(Loadings_450k_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=2)+theme_bw()+fillscale_sampsite

Loadings_EPIC$array.id<-rownames(Loadings_EPIC)
Loadings_EPIC_meta<-merge(Loadings_EPIC, sampleinfo_EPIC, by="array.id")
epic<-ggplot(Loadings_EPIC_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=2)+theme_bw()+fillscale_sampsite
grid.arrange(k450, epic)

k450<-ggplot(Loadings_450k_meta, aes(PC1, PC2, fill=diagnosis_grouped))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis
epic<-ggplot(Loadings_EPIC_meta, aes(PC1, PC2, fill=diagnosis_grouped))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis
grid.arrange(k450, epic)




#' ### PCA split by tissues
sampleinfo_TI<-sampleinfo[which(sampleinfo$sample.site=="TI"),]
ibd_combo_TI<-ibd_combo[,which(colnames(ibd_combo)%in%sampleinfo_TI$array.id)]
identical(colnames(ibd_combo_TI),sampleinfo_TI$array.id)

sampleinfo_C<-sampleinfo[which(sampleinfo$sample.site!="TI"),]
ibd_combo_C<-ibd_combo[,which(colnames(ibd_combo)%in%sampleinfo_C$array.id)]
identical(colnames(ibd_combo_C),sampleinfo_C$array.id)

pca_res_TI <- prcomp(t(ibd_combo_TI))
Loadings_TI<-as.data.frame(pca_res_TI$x)

pca_res_C <- prcomp(t(ibd_combo_C))
Loadings_C<-as.data.frame(pca_res_C$x)

## PC vs PC plot
Loadings_TI$array.id<-rownames(Loadings_TI)
Loadings_TI_meta<-merge(Loadings_TI, sampleinfo_TI, by="array.id")
TI<-ggplot(Loadings_TI_meta, aes(PC1, PC2, fill=diagnosis_grouped))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis

Loadings_C$array.id<-rownames(Loadings_C)
Loadings_C_meta<-merge(Loadings_C, sampleinfo_C, by="array.id")
C<-ggplot(Loadings_C_meta, aes(PC1, PC2, fill=diagnosis_grouped))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis
grid.arrange(TI, C)


## taking the fit of PC1 and PC2 and seeing if it assocaites to anything
colon<-Loadings_all_meta[which(Loadings_all_meta$sample.site_grouped=="Colon"),]

z<-lm(colon$PC1~colon$PC2)
intercept=z$coefficients[1]
slope=z$coefficients[2]

colon$pc_trend<-sapply(1:nrow(colon), function(x){
  (slope*colon$PC1[x]+intercept)
})



ileum<-Loadings_all_meta[which(Loadings_all_meta$sample.site_grouped=="Small\nIntestine"),]
z<-lm(ileum$PC1~ileum$PC2)
intercept=z$coefficients[1]
slope=z$coefficients[2]

ileum$pc_trend<-sapply(1:nrow(ileum), function(x){
  (slope*ileum$PC1[x]+intercept)
})


ggplot(Loadings_all_meta, aes(PC1, PC2, fill=sample.site_grouped))+geom_point(shape=21,size=3)+theme_bw()+
  th+scale_fill_manual(values=c("#a6d96a","cornflowerblue"), name="Sample Site")+stat_smooth(method="lm", color="black")

ggsave(here("DNAm/figs","PC1_PC2_linear_fit.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_linear_fit.jpeg"), width = 7.5, height = 6)


ggplot(colon, aes(PC1, PC2, fill=pc_trend))+geom_point(shape=21,size=3)+theme_bw()
ggplot(ileum, aes(PC1, PC2, fill=pc_trend))+geom_point(shape=21,size=3)+theme_bw()


summary(aov(colon$pc_trend ~ colon$inflammation+colon$diagnosis_grouped))
summary(aov(ileum$pc_trend ~ ileum$inflammation+ileum$diagnosis_grouped))

summary(aov(colon$pc_trend ~ colon$diagnosis_grouped+colon$inflammation))
summary(aov(ileum$pc_trend ~ ileum$diagnosis_grouped+ileum$inflammation))

summary(aov(colon$pc_trend ~ colon$array.type))
summary(aov(ileum$pc_trend ~ ileum$array.type))


cor.test(ileum$pc_trend , ileum$age)
cor.test(colon$pc_trend , colon$age)


summary(aov(colon$pc_trend ~ colon$sex))
summary(aov(ileum$pc_trend ~ ileum$sex))

ggplot(colon, aes(inflammation,pc_trend, fill=inflammation))+
  geom_violin(alpha=0.5)+geom_boxplot()+th+theme_bw()+scale_fill_manual(values=c("lightgrey","#e34a33"))+facet_wrap(~diagnosis_grouped)

ggplot(ileum, aes(inflammation,pc_trend, fill=inflammation))+
  geom_violin(alpha=0.5)+geom_boxplot()+th+theme_bw()+scale_fill_manual(values=c("lightgrey","#e34a33"))+facet_wrap(~diagnosis_grouped)


ggplot(colon, aes(age,pc_trend, fill=age))+
  geom_point()+th+theme_bw()+stat_smooth(method="lm", color="black")

ggplot(ileum, aes(age,pc_trend, fill=age))+
  geom_point()+th+theme_bw()+stat_smooth(method="lm", color="black")


ggplot(colon, aes(diagnosis_grouped,pc_trend, fill=diagnosis_grouped))+
  geom_violin(alpha=0.5)+geom_boxplot()+th+theme_bw()+fillscale_diagnosis

ggplot(ileum, aes(diagnosis_grouped,pc_trend, fill=diagnosis_grouped))+
  geom_violin(alpha=0.5)+geom_boxplot()+th+theme_bw()+fillscale_diagnosis

inf<-ggplot(Loadings_all_meta, aes(PC1, PC2, fill=inflammation))+geom_point(shape=21,size=2)+theme_bw()+
  scale_fill_manual(values=c("grey", "red"))+th+theme(plot.margin=unit(c(0.25,0.15,0.25,0.175),"cm"))
age<-ggplot(Loadings_all_meta, aes(PC1, PC2, fill=age))+geom_point(shape=21,size=2)+theme_bw()+
  th+theme(plot.margin=unit(c(0.25,1.35,0.25,0.2),"cm"))
sex<-ggplot(Loadings_all_meta, aes(PC1, PC2, fill=sex))+geom_point(shape=21,size=2)+theme_bw()+
  scale_fill_manual(values=c("#e8702a", "#bcd2d0"))+th+theme(plot.margin=unit(c(0.25,1.55,0.25,0.25),"cm"))
array<-ggplot(Loadings_all_meta, aes(PC1, PC2, fill=array.type))+geom_point(shape=21,size=2)+theme_bw()+
  scale_fill_manual(values=c("darkblue","springgreen4"), name="Array Type")+th+theme(plot.margin=unit(c(0.25,0.4,0.25,0.25),"cm"))

grid.arrange(inf, age, array, sex)

ggsave(here("DNAm/figs","PC1_PC2_spreadexplained.pdf"),grid.arrange(inf, age, array, sex), width = 9, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2__spreadexplained.jpeg"),grid.arrange(inf, age, array, sex), width = 9, height = 6)


################################# PCA down sampled to one sample per person
pltpltplt<-lapply(1:16, function(sed){
    sampleinfo_one_only<-lapply(1:length(unique(sampleinfo$case.no)), function(x){
      case<-sampleinfo[which(sampleinfo$case.no==unique(sampleinfo$case.no)[x]),]
      set.seed(sed)
      case[sample(1:nrow(case), 1),]
    })
    
    sampleinfo_one_only<-do.call(rbind, sampleinfo_one_only)
    
    ibd_combo_one_only<-ibd_combo[,which(colnames(ibd_combo)%in%sampleinfo_one_only$array.id)]
    identical(colnames(ibd_combo_one_only),sampleinfo_one_only$array.id)
    
    pca_res_onlyone <- prcomp(t(ibd_combo_one_only))
    Loadings_onlyone<-as.data.frame(pca_res_onlyone$x)
    
    ## PC vs PC plot
    Loadings_onlyone$array.id<-rownames(Loadings_onlyone)
    Loadings_onlyone_meta<-merge(Loadings_onlyone, sampleinfo_one_only, by="array.id")
    ggplot(Loadings_onlyone_meta, aes(PC1, PC2, fill=sample.site))+geom_point(shape=21,size=2)+theme_bw()+fillscale_sampsite})

do.call("grid.arrange", c(pltpltplt, ncol=4))



########## so it isn't consistent in people across tissues btu maybe with rescopes in the same tissue it is consistent
sampleinfo$diagnosis_grouped<-as.factor(sampleinfo$diagnosis)
levels(sampleinfo$diagnosis_grouped)<-c("CD","Control","CD","UC","UC")
sampleinfo$diagnosis_grouped<-factor(sampleinfo$diagnosis_grouped, levels=c("Control", "CD", "UC"))


pca_res_all <- prcomp(t(ibd_combo))
Loadings_all<-as.data.frame(pca_res_all$x)

Loadings_all$array.id<-rownames(Loadings_all)
Loadings_all_meta<-merge(Loadings_all, sampleinfo, by="array.id")
ggplot(Loadings_all_meta, aes(PC1, PC2, fill=diagnosis_grouped))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis

ggplot(Loadings_all_meta, aes(PC1, PC2, fill=sampling.time.point))+geom_point(shape=21,size=2, alpha=0.2)+
  geom_line(aes(PC1,PC2, group=sample_ID))+
  theme_bw()+scale_fill_manual(values=c("grey93","blue"))

ggsave(here("DNAm/figs","PC1_PC2_distance_rescopes.pdf"), width = 10, height = 8)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_distance_rescopes.jpeg"), width = 10, height = 8)






#'## R Session Info
sessionInfo()



