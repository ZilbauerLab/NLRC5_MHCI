#'---
#'title: Differential DNAm 
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(limma)
library(here)

options(stringsAsFactors = FALSE)


#' functions
source(here("general_functions","00_pretty_plots.R"))
source(here("general_functions","00_Heat_scree_plot_generic.R"))
source(here("general_functions","00_DNAm_volcano.R"))



generic_boxplot<-function(cpgs){
  beta<-combat_organoid_Beta[cpgs,]
  beta<-melt(beta)
  if(length(cpgs)==1){
    beta$Var2<-rownames(beta)
    beta$Var1<-cpgs}
  plt<-merge(beta, epic.organoid_exclude_other, by.x="Var2", by.y="array.id")
  ggplot(plt, aes(diagnosis, value, fill=diagnosis))+
    geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.15, shape=21, size=1,color="black")+
    facet_grid(Var1~sample.site)+fillscale_diagnosis+th+theme_bw()+ylim(0,1)+ylab("DNA Methylation")+xlab("Diagnosis")}



#' ## Load DNAm data
load(file=paste(here("DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))
          #combat_organoid_Beta<-combat_organoid_Beta[1:1000,]

## fix 202 and 203
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605870061_R03C01")]<-"hold"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605880096_R04C01")]<-"205605870061_R03C01"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="hold")]<-"205605880096_R04C01"



#' ## convert to M values for stats
Mval<-function(beta) log2(beta/(1-beta))
ibd_organoid_Mval = apply(combat_organoid_Beta, 1, Mval) # need mvalues for combat
ibd_organoid_Mval = as.data.frame(ibd_organoid_Mval)
ibd_organoid_Mval = t(ibd_organoid_Mval)

# # only plot a subset fo CpGs for speed
# set.seed(1)
# mval_melted<- melt(ibd_organoid_Mval[sample(nrow(ibd_organoid_Mval),10000),])
# mval_Plot<-mval_melted[which(!(is.na(mval_melted$value))),]
# colnames(mval_Plot)<-c("CpG","ID","Mval")
# mval_Plot<-merge(mval_Plot,epic.organoid, by.x="ID", by.y="array.id")
# 
# ggplot(mval_Plot, aes(Mval, group=as.character(sample_ID), color=as.character(sample.site)))+
#   geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm M Value")
# 
# ggsave(here("DNAm/figs","Mvalue_organoid_distribution_combined.pdf"), w=10, h=5)
# ggsave(here("DNAm/figs/jpeg","Mvalue_organoid_distribution_combined.jpeg"),  w=10, h=5)
# 



### load update diagnosis
epic.samples<-read.table(here("DNAm/data","AllEPICOrganoid_MHCIannotated_UpdatedSampleInfo_31Oct22.txt"), header=T, sep="\t")

head(epic.organoid_combined[,c("case.no","array.id" , "sample.site" ,"array.type","sentrix.pos" , 
                               "sample_ID","sentrix_ID","det_pval","passage.or.rescope.no_numeric", "Sample_Name",
                               "Wnt.type","condition","Biobank.Rachel.Replicates","batch")])

head(epic.samples[,c("array.id","diagnosis", "sex", "age","sample.id","numeric.passage")])

epic.organoid_combined_update<-merge(epic.organoid_combined[,c("array.id" , "sample.site" ,"array.type","sentrix.pos" , 
                                                               "sample_ID","sentrix_ID","det_pval","passage.or.rescope.no_numeric", "Sample_Name",
                                                               "Wnt.type","condition","Biobank.Rachel.Replicates","batch")],
                                     epic.samples[,c("case.no", "array.id","diagnosis", "sex", "age","sample.id","numeric.passage",
                                                     "Disease.Severity","DUO.Inflammation", "TI.Inflammation", "SC.Inflammation", 
                                                     "Biologics", "Surgery","AZA", "Perianal.Disease")],
                                     by="array.id")
dim(epic.organoid_combined_update)
epic.organoid_combined_update[grep("202|203", epic.organoid_combined_update$case.no),]


## for reporting purposes 6 samples used are rescopes, the rest are diagnosis
table(epic.samples[grep("RE|re",epic.samples$RE.sample.name),]$tissue)
#Duo  SC  TI
# 2   1   3



epic.organoid_exclude_other<-epic.organoid_combined_update
epic.organoid_exclude_other$sample.site<-as.factor(epic.organoid_exclude_other$sample.site)
levels(epic.organoid_exclude_other$sample.site)<-c("DUO","DUO","SC","TI")

## this is all in seperate code now
#' 
#' 
#' ## check over rescopes
#' rescopes<-epic.organoid_combined_update[grep("RE",epic.organoid_combined_update$passage.or.rescope.no),]
#' epic.organoid_combined_update[which(epic.organoid_combined_update$sample.id%in%rescopes$sample.id),]
#' # T091 SC and T116 SC I have an original and rescope
#' 
#' #' ## Exclude other.GI
#' epic.organoid_exclude_other<-epic.organoid_combined_update[which(!(epic.organoid_combined_update$diagnosis%in%c("Other.GI"))),]
#' epic.organoid_exclude_other$diagnosis<-as.factor(epic.organoid_exclude_other$diagnosis)
#' levels(epic.organoid_exclude_other$diagnosis)
#' 
#' #'## remove IFNg, TNFa, LPS, AZA treated as well as differentiated but left in frozen biopsies
#' epic.organoid_exclude_other<-epic.organoid_exclude_other[which(!(epic.organoid_exclude_other$condition%in%c("D","IFNg","TNFa","LPS","AZA-3W","AZA-72H" ))),]
#' 
#' ## check over rescopes
#' rescopes<-epic.organoid_exclude_other[grep("RE",epic.organoid_exclude_other$passage.or.rescope.no),]
#' epic.organoid_exclude_other[which(epic.organoid_exclude_other$sample.id%in%rescopes$sample.id),]
#' # T091 SC and T116 SC I have an original and rescope
#' epic.organoid_exclude_other[grep("T091|T116",epic.organoid_exclude_other$case.no),]
#' # will prioritize keeping the diagnosis sample for T091 over the lowest passage (RE P2 original P4) for T116 original is the lowest passage
#' epic.organoid_exclude_other<-epic.organoid_exclude_other[-which(epic.organoid_exclude_other$Sample_Name=="T091 SC RE1.P2"),]
#' 
#' #' ## use only sample with fewest passages
#' #' when there are mulitple passages from the same individual and sample site
#' samples<-unique(epic.organoid_exclude_other$sample.id)
#' epic.organoid_exclude_other<-lapply(1:length(samples), function(x){
#'   samp<-epic.organoid_exclude_other[which(epic.organoid_exclude_other$sample.id==samples[x]),]
#'   samp[which(samp$passage.or.rescope.no_numeric==min(samp$passage.or.rescope.no_numeric)),]
#' })
#' epic.organoid_exclude_other<-do.call(rbind, epic.organoid_exclude_other)
#' 
#' #'## remove cross batch replicates
#' dupli<-epic.organoid_exclude_other[duplicated(epic.organoid_exclude_other$sample.id),]
#' epic.organoid_exclude_other[which(epic.organoid_exclude_other$sample.id%in%dupli$sample.id),]
#' 
#' #'#### remove the replicates where one is  frozen biopsy
#' epic.organoid_exclude_other<-epic.organoid_exclude_other[which(!(epic.organoid_exclude_other$Sample_Name%in%c("T074 RE3 DUO P2 FB","T421 DUO P2 FB"))),]
#' #'####otherwise keep the new
#' dupli<-epic.organoid_exclude_other[duplicated(epic.organoid_exclude_other$sample.id),]
#' epic.organoid_exclude_other[which(epic.organoid_exclude_other$sample.id%in%dupli$sample.id),]
#' rm_old_dup<-epic.organoid_exclude_other[which(epic.organoid_exclude_other$sample.id%in%dupli$sample.id),]
#' rm_old_dup<-rm_old_dup$array.id[which(rm_old_dup$batch=="old")]
#' epic.organoid_exclude_other<-epic.organoid_exclude_other[which(!(epic.organoid_exclude_other$array.id%in%rm_old_dup)),]
#' dupli<-epic.organoid_exclude_other[duplicated(epic.organoid_exclude_other$sample.id),]
#' rm_old_dup<-epic.organoid_exclude_other[which(epic.organoid_exclude_other$sample.id%in%dupli$sample.id),]
#' rm_old_dup<-rm_old_dup$array.id[which(rm_old_dup$batch=="new")]
#' epic.organoid_exclude_other<-epic.organoid_exclude_other[which(!(epic.organoid_exclude_other$array.id%in%rm_old_dup)),]
#' 
#' 
#' #'#### no duplicates remain
#' dupli<-epic.organoid_exclude_other[duplicated(epic.organoid_exclude_other$sample.id),]
#' epic.organoid_exclude_other[which(epic.organoid_exclude_other$sample.id%in%dupli$sample.id),]
#' nrow(epic.organoid_exclude_other) #316
#' length(unique(epic.organoid_exclude_other$sample.id))
#' table(epic.organoid_exclude_other$sample.site)
#' 
#' ## Remove fetal and neonatal
#' table(epic.organoid_exclude_other$Biobank.Rachel.Replicates)
#' epic.organoid_exclude_other<-epic.organoid_exclude_other[which(epic.organoid_exclude_other$Biobank.Rachel.Replicates!="Natasha"),]
#' 
#' 
# 
#' #' ## some updates in diagnosis from Matt
#' epic.organoid_exclude_other[which(epic.organoid_exclude_other$case.no=="T472"),]$diagnosis<-"CD"
#' epic.organoid_exclude_other[which(epic.organoid_exclude_other$case.no=="T027"),]$diagnosis<-"CD"
#' epic.organoid_exclude_other[which(epic.organoid_exclude_other$case.no=="T490"),]$TI.Inflammation<-"Macroscopic"
#' epic.organoid_exclude_other[which(epic.organoid_exclude_other$case.no=="T490"),]$diagnosis<-"CD"
#' 








ibd_organoid_Mval<-ibd_organoid_Mval[,which(colnames(ibd_organoid_Mval)%in%epic.organoid_exclude_other$array.id)]
ibd_organoid_Mval<-ibd_organoid_Mval[,match(epic.organoid_exclude_other$array.id,colnames(ibd_organoid_Mval) )]
identical(colnames(ibd_organoid_Mval),epic.organoid_exclude_other$array.id)

combat_organoid_Beta<-combat_organoid_Beta[,which(colnames(combat_organoid_Beta)%in%epic.organoid_exclude_other$array.id)]
combat_organoid_Beta<-combat_organoid_Beta[,match(epic.organoid_exclude_other$array.id,colnames(combat_organoid_Beta) )]
identical(colnames(combat_organoid_Beta),epic.organoid_exclude_other$array.id)

## sample size for comparisons
tapply(epic.organoid_exclude_other$array.id, list(epic.organoid_exclude_other$sample.site, epic.organoid_exclude_other$diagnosis), function(x) length(unique(x)))


#' ## diagnosis simple
epic.organoid_exclude_other$diagnosis_simple<-as.factor(epic.organoid_exclude_other$diagnosis)
levels(epic.organoid_exclude_other$diagnosis_simple)<-c("IBD","Control","IBD","IBD")

#'### Plot of passages used
ggplot(epic.organoid_exclude_other, aes(as.factor(passage.or.rescope.no_numeric), fill=as.factor(passage.or.rescope.no_numeric)))+geom_bar(color="black")+
  theme_bw()+theme(legend.position = "none", 
                   legend.title=element_text(size=10),
                   legend.text=element_text(size=8))+xlab("Passage Number")+
  scale_fill_manual(values=pass_col,name="Passage\nNumber")+ylab("Sample Number")+th
ggsave(here("DNAm/figs","organoid_passage_num_EWAS_ThreeBatch.pdf"), width = 6, height = 5)
ggsave(here("DNAm/figs/jpeg","organoid_passage_num_EWAS_ThreeBatch.jpeg"), width = 6, height = 5)



#' ## subset tissues
sampleinfo_organoid_exclude_TI<-epic.organoid_exclude_other[which(epic.organoid_exclude_other$sample.site=="TI"),]
ibd_Mval_organoid_exclude_TI<-ibd_organoid_Mval[,which(colnames(ibd_organoid_Mval)%in%sampleinfo_organoid_exclude_TI$array.id)]
ibd_combo_organoid_exclude_TI<-combat_organoid_Beta[,which(colnames(combat_organoid_Beta)%in%sampleinfo_organoid_exclude_TI$array.id)]
identical(colnames(ibd_Mval_organoid_exclude_TI),sampleinfo_organoid_exclude_TI$array.id)

sampleinfo_organoid_exclude_SC<-epic.organoid_exclude_other[which(epic.organoid_exclude_other$sample.site=="SC"),]
ibd_Mval_organoid_exclude_SC<-ibd_organoid_Mval[,which(colnames(ibd_organoid_Mval)%in%sampleinfo_organoid_exclude_SC$array.id)]
ibd_combo_organoid_exclude_SC<-combat_organoid_Beta[,which(colnames(combat_organoid_Beta)%in%sampleinfo_organoid_exclude_SC$array.id)]
identical(colnames(ibd_Mval_organoid_exclude_SC),sampleinfo_organoid_exclude_SC$array.id)
identical(colnames(ibd_combo_organoid_exclude_SC),sampleinfo_organoid_exclude_SC$array.id)

sampleinfo_organoid_exclude_duo<-epic.organoid_exclude_other[which(epic.organoid_exclude_other$sample.site=="DUO"),]
ibd_Mval_organoid_exclude_duo<-ibd_organoid_Mval[,which(colnames(ibd_organoid_Mval)%in%sampleinfo_organoid_exclude_duo$array.id)]
ibd_combo_organoid_exclude_duo<-combat_organoid_Beta[,which(colnames(combat_organoid_Beta)%in%sampleinfo_organoid_exclude_duo$array.id)]
identical(colnames(ibd_Mval_organoid_exclude_duo),sampleinfo_organoid_exclude_duo$array.id)
identical(colnames(ibd_combo_organoid_exclude_duo),sampleinfo_organoid_exclude_duo$array.id)



#' ## meta correlation between potential covariates
meta_categorical <- epic.organoid_exclude_other[, c("diagnosis","sample.site","sex")]  # input column numbers in meta that contain categorical variables
meta_continuous <- as.data.frame(epic.organoid_exclude_other[, c("age", "passage.or.rescope.no_numeric")] ) # input column numbers in meta that contain continuous variables
colnames(meta_categorical) <- c("Diagnosis", "Sample Site","Sex")
colnames(meta_continuous) <- c("Age", "Passage")

Meta_correlation<-lapply(1:ncol(meta_continuous), function(x) as.numeric(meta_continuous[,x]))
Meta_correlation<-as.data.frame(do.call(cbind, Meta_correlation))
colnames(Meta_correlation)<-colnames(meta_continuous)

cor_meta<-lapply(1:ncol(meta_continuous), function(con1) {
  sapply(1:ncol(meta_continuous), function(con2){
    cor.test(meta_continuous[,con1] , meta_continuous[,con2])$p.value})
})
cor_meta<-as.data.frame(do.call(rbind, cor_meta))
colnames(cor_meta)<-colnames(meta_continuous)

aov_meta<-lapply(1:ncol(meta_categorical), function(cat) {
  sapply(1:ncol(meta_continuous), function(con1){
    summary(aov(meta_continuous[, con1] ~ meta_categorical[, cat]))[[1]]$"Pr(>F)"[1] })
})
aov_meta<-as.data.frame(do.call(rbind, aov_meta))
colnames(aov_meta)<-colnames(meta_continuous)

suppressWarnings(chi_meta<-lapply(1:ncol(meta_categorical), function(cat1) {
  sapply(1:ncol(meta_categorical), function(cat2){
    chisq.test(table(meta_categorical[,cat1], meta_categorical[,cat2]))$p.value})
}))

chi_meta<-as.data.frame(do.call(rbind, chi_meta))
colnames(chi_meta)<-colnames(meta_categorical)

#add continious
chi_meta<-cbind(chi_meta,aov_meta)
con_t<-t(aov_meta)
colnames(con_t)<-colnames(meta_categorical)
chi_meta<-rbind(chi_meta,cbind(con_t,cor_meta))

chi_meta$id<-c(colnames(meta_categorical),"Age","Passage")
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
  fillscale_pval

ggsave(here("DNAm/figs","organoid_meta_correlation_ThreeBatch.pdf"), width = 6, height = 5)
ggsave(here("DNAm/figs/jpeg","organoid_correlation_ThreeBatch.jpeg"), width = 6, height = 5)

#' ### proportion of diagnosis and sex and age
pro_sex<-as.data.frame.matrix(table(epic.organoid_exclude_other$diagnosis, epic.organoid_exclude_other$sex))
pro_sex$diagnosis<-rownames(pro_sex)
pro_sex
pro_sex<-melt(pro_sex)

ggplot(pro_sex, aes(diagnosis,value, fill=variable))+geom_bar(stat="identity", position="dodge")+theme_bw()

pro_age<-as.data.frame.matrix(table(epic.organoid_exclude_other$diagnosis, epic.organoid_exclude_other$age))
pro_age$diagnosis<-rownames(pro_age)
pro_age
pro_age<-melt(pro_age)

ggplot(pro_age, aes(diagnosis,value))+geom_bar(stat="identity", position="dodge")+theme_bw()



#######################
#'##  IBD vs controls 
#' only relevant in SC
#######################
DMA<-function(mval,beta,meta){
  print(paste("Testing in ", ncol(beta), " individuals and only ", unique(meta$sample.site), " samples", sep=""))
  # p value
  mod<-model.matrix(~as.factor(diagnosis_simple)+as.factor(sex)+age, data=meta)
  fit <- lmFit(mval, mod)
  ebfit <- eBayes(fit)
  
  # covariate adjusted beta values
  avebeta.lm<-apply(beta, 1, function(x){
    lm(x~as.factor(sex)+age, data=meta)})
  residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
  colnames(residuals)<-colnames(beta)
  adj.residuals<-residuals+matrix(apply(beta, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
  
  delta_beta<-sapply(1:nrow(adj.residuals), function(x) {
    group_means<-tapply(adj.residuals[x,], meta$diagnosis_simple, mean)
    group_means[1]-group_means[2] # IBD-CTRL
  })
  
  limmares <- data.frame(p.value=ebfit$p.value[,"as.factor(diagnosis_simple)Control"], delta_beta=delta_beta)
  limmares}

diff_meth_TI_organoid<-DMA(ibd_Mval_organoid_exclude_TI,ibd_combo_organoid_exclude_TI, sampleinfo_organoid_exclude_TI)
diff_meth_SC_organoid<-DMA(ibd_Mval_organoid_exclude_SC,ibd_combo_organoid_exclude_SC, sampleinfo_organoid_exclude_SC)
diff_meth_duo_organoid<-DMA(ibd_Mval_organoid_exclude_duo,ibd_combo_organoid_exclude_duo, sampleinfo_organoid_exclude_duo)

diff_meth_TI_organoid$adjusted_p<-p.adjust(diff_meth_TI_organoid$p.value, method = "fdr", n = nrow(diff_meth_TI_organoid))
diff_meth_SC_organoid$adjusted_p<-p.adjust(diff_meth_SC_organoid$p.value, method = "fdr", n = nrow(diff_meth_SC_organoid))
diff_meth_duo_organoid$adjusted_p<-p.adjust(diff_meth_duo_organoid$p.value, method = "fdr", n = nrow(diff_meth_duo_organoid))

save(diff_meth_TI_organoid,diff_meth_SC_organoid,diff_meth_duo_organoid, file=here("DNAm/data","organoid_diff_DNAm_simple_casectrl_ThreeBatch.RData"))


#######################
#' ## disease subset specific grouping
#######################
DMA_grouped<-function(mval,beta,meta){
  #remove IBD-U
  meta<-meta[which(meta$diagnosis!="IBD-U"),]
  beta<-beta[,which(colnames(beta)%in%meta$array.id)]
  mval<-mval[,which(colnames(mval)%in%meta$array.id)]
  print(paste("Testing in ", ncol(beta), " individuals and only ", unique(meta$sample.site), " samples", sep=""))
  
  # p value
  meta$diagnosis<-as.factor(as.character(meta$diagnosis))
  meta$sex<-as.factor(meta$sex)
  
  mod<-model.matrix(~0+diagnosis+sex+age, data=meta)
  fit <- lmFit(mval, mod)
  #ebfit <- eBayes(fit)
  
  # construct the contrast matrix
  contrastMatrix <- makeContrasts(
    ctrl_CD = diagnosisControl - diagnosisCD,
    ctrl_UC = diagnosisControl - diagnosisUC,
    CD_UC = diagnosisCD - diagnosisUC,
    levels = mod
  )
  
  contrastFit<-contrasts.fit(fit, contrastMatrix)
  contrastFitEb <- eBayes(contrastFit) 
  contrastCpG_CD <- topTable(contrastFitEb, coef="ctrl_CD",number=Inf)
  contrastCpG_UC <- topTable(contrastFitEb, coef="ctrl_UC",number=Inf)
  contrastCpG_CDUC <- topTable(contrastFitEb, coef="CD_UC",number=Inf)
  
  # covariate adjusted beta values
  avebeta.lm<-apply(beta, 1, function(x){
    lm(x~as.factor(sex)+age, data=meta)})
  residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
  colnames(residuals)<-colnames(beta)
  adj.residuals<-residuals+matrix(apply(beta, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
  
  delta_beta<-lapply(1:nrow(adj.residuals), function(x) {
    group_means<-tapply(adj.residuals[x,], meta$diagnosis, mean)
    c(group_means[2]-group_means[1],group_means[2]-group_means[3], group_means[3]-group_means[1])
  })
  
  limmares<-as.data.frame(do.call(rbind, delta_beta))
  colnames(limmares)<-c("db_ctrl_CD","db_ctrl_UC","db_UC_CD")
  limmares$CpG<-rownames(adj.residuals)
  
  
  p.value_UC<-data.frame(p.value_UC=contrastCpG_UC$P.Value, CpG=rownames(contrastCpG_UC))
  p.value_CD<-data.frame(p.value_CD=contrastCpG_CD$P.Value, CpG=rownames(contrastCpG_CD))
  p.value_CDUC<-data.frame(p.value_CDUC=contrastCpG_CDUC$P.Value, CpG=rownames(contrastCpG_CDUC))
  
  
  pval<-merge(p.value_CDUC,p.value_UC, by="CpG")
  pval<-merge(pval,p.value_CD, by="CpG")
  limmares<-merge(limmares,pval, by="CpG")
  
  limmares$adjusted_p_CDUC<-p.adjust(limmares$p.value_CDUC, method = "fdr", n = nrow(limmares))
  limmares$adjusted_p_UC<-p.adjust(limmares$p.value_UC, method = "fdr", n = nrow(limmares))
  limmares$adjusted_p_CD<-p.adjust(limmares$p.value_CD, method = "fdr", n = nrow(limmares))
  
  limmares}

diff_meth_TI_grouped<-DMA_grouped(ibd_Mval_organoid_exclude_TI,ibd_combo_organoid_exclude_TI, sampleinfo_organoid_exclude_TI)
diff_meth_SC_grouped<-DMA_grouped(ibd_Mval_organoid_exclude_SC,ibd_combo_organoid_exclude_SC, sampleinfo_organoid_exclude_SC)
diff_meth_duo_grouped<-DMA_grouped(ibd_Mval_organoid_exclude_duo,ibd_combo_organoid_exclude_duo, sampleinfo_organoid_exclude_duo)

save(diff_meth_TI_grouped, diff_meth_SC_grouped,diff_meth_duo_grouped, file=here("DNAm/data/organoid_diff_DNAm_grouped_casectrl_ThreeBatch.RData"))




####################
#' # Grouping UC with controls in TI makes biological sense
####################

# IBD vs Ctrl (done already)
# UC vs Ctrl in SC (done already)
# CD vs Ctrl in SC (done already)
# CD vs UC+Ctrl in TI (done below)

mval<-ibd_Mval_organoid_exclude_TI
beta<-ibd_combo_organoid_exclude_TI
meta<-sampleinfo_organoid_exclude_TI
#remove IBD-U
meta<-meta[which(meta$diagnosis!="IBD-U"),]
beta<-beta[,which(colnames(beta)%in%meta$array.id)]
mval<-mval[,which(colnames(mval)%in%meta$array.id)]
print(paste("Testing in ", ncol(beta), " individuals and only ", unique(meta$sample.site), " samples", sep=""))
# p value
meta$diagnosis<-as.factor(as.character(meta$diagnosis))
meta$diagnosis_TI_split<-meta$diagnosis
levels(meta$diagnosis_TI_split)<-c("CD","Control_UC","Control_UC")

meta$sex<-as.factor(meta$sex)

mod<-model.matrix(~diagnosis_TI_split+sex+age, data=meta)
fit <- lmFit(mval, mod)
ebfit <- eBayes(fit)

# covariate adjusted beta values
avebeta.lm<-apply(beta, 1, function(x){
  lm(x~as.factor(sex)+age, data=meta)})
residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
colnames(residuals)<-colnames(beta)
adj.residuals<-residuals+matrix(apply(beta, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))

delta_beta<-sapply(1:nrow(adj.residuals), function(x) {
  group_means<-tapply(adj.residuals[x,], meta$diagnosis_TI_split, mean)
  group_means[2]-group_means[1] # Control_UC - CD
})

diff_meth_TI_UCCTRL_organoid <- data.frame(p.value=ebfit$p.value[,"diagnosis_TI_splitControl_UC"], delta_beta=delta_beta)

save(diff_meth_TI_UCCTRL_organoid, file=here("DNAm/data/diff_DNAm_TI_UCCTRL_organoid_ThreeBatch.RData"))



#### DUO same model as TI
mval<-ibd_Mval_organoid_exclude_duo
beta<-ibd_combo_organoid_exclude_duo
meta<-sampleinfo_organoid_exclude_duo
#remove IBD-U
meta<-meta[which(meta$diagnosis!="IBD-U"),]
beta<-beta[,which(colnames(beta)%in%meta$array.id)]
mval<-mval[,which(colnames(mval)%in%meta$array.id)]
print(paste("Testing in ", ncol(beta), " individuals and only ", unique(meta$sample.site), " samples", sep=""))
# p value
meta$diagnosis<-as.factor(as.character(meta$diagnosis))
meta$diagnosis_DUO_split<-meta$diagnosis
levels(meta$diagnosis_DUO_split)<-c("CD","Control_UC","Control_UC")

meta$sex<-as.factor(meta$sex)

mod<-model.matrix(~diagnosis_DUO_split+sex+age, data=meta)
fit <- lmFit(mval, mod)
ebfit <- eBayes(fit)

# covariate adjusted beta values
avebeta.lm<-apply(beta, 1, function(x){
  lm(x~as.factor(sex)+age, data=meta)})
residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
colnames(residuals)<-colnames(beta)
adj.residuals<-residuals+matrix(apply(beta, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))

delta_beta<-sapply(1:nrow(adj.residuals), function(x) {
  group_means<-tapply(adj.residuals[x,], meta$diagnosis_DUO_split, mean)
  group_means[2]-group_means[1] # Control_UC - CD
})

diff_meth_DUO_UCCTRL_organoid <- data.frame(p.value=ebfit$p.value[,"diagnosis_DUO_splitControl_UC"], delta_beta=delta_beta)

save(diff_meth_DUO_UCCTRL_organoid, file=here("DNAm/data/diff_DNAm_DUO_UCCTRL_organoid_ThreeBatch.RData"))



#' ### Load results

# IBD vs Ctrl
load(here("DNAm/data/organoid_diff_DNAm_simple_casectrl_ThreeBatch.RData"))
# UC vs Ctrl in SC
# CD vs Ctrl in SC
load(here("DNAm/data/organoid_diff_DNAm_grouped_casectrl_ThreeBatch.RData"))
# CD vs UC+Ctrl in TI
load(here("DNAm/data/diff_DNAm_TI_UCCTRL_organoid_ThreeBatch.RData"))
diff_meth_TI_UCCTRL_organoid$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL_organoid$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL_organoid))
# CD vs UC+Ctrl in DUO
load(file=here("DNAm/data/diff_DNAm_DUO_UCCTRL_organoid_ThreeBatch.RData"))
diff_meth_DUO_UCCTRL_organoid$adjusted_p<-p.adjust(diff_meth_DUO_UCCTRL_organoid$p.value, method = "fdr", n = nrow(diff_meth_DUO_UCCTRL_organoid))



#' ## Significant CpGs numbers
dB=0.05
fdr=0.05

# SC IBD
stat_hits<-as.data.frame(diff_meth_SC_organoid)[which(diff_meth_SC_organoid$adjusted_p<=fdr),]
SC_IBD<-stat_hits[which(abs(stat_hits$delta_beta)>=dB),]
SC_IBD$CpG<-rownames(SC_IBD)
nrow(SC_IBD)

# SC CD
stat_hits<-as.data.frame(diff_meth_SC_grouped)[which(diff_meth_SC_grouped$adjusted_p_CD<=fdr),]
SC_CD<-stat_hits[which(abs(stat_hits$db_ctrl_CD)>=dB),]
nrow(SC_CD)

# SC UC
stat_hits<-as.data.frame(diff_meth_SC_grouped)[which(diff_meth_SC_grouped$adjusted_p_UC<=fdr),]
SC_UC<-stat_hits[which(abs(stat_hits$db_ctrl_UC)>=dB),]
nrow(SC_UC)

# TI CD
stat_hits<-as.data.frame(diff_meth_TI_UCCTRL_organoid)[which(diff_meth_TI_UCCTRL_organoid$adjusted_p<=fdr),]
TI_sig_UCCTRL<-stat_hits[which(abs(stat_hits$delta_beta)>=dB),]
TI_sig_UCCTRL$CpG<-rownames(TI_sig_UCCTRL)
nrow(TI_sig_UCCTRL)

          
          # DUO CD
          stat_hits<-as.data.frame(diff_meth_duo_grouped)[which(diff_meth_duo_grouped$adjusted_p_CD<=fdr),]
          DUO_CD<-stat_hits[which(abs(stat_hits$db_ctrl_CD)>=dB),]
          nrow(DUO_CD)
          
          # DUO UC
          stat_hits<-as.data.frame(diff_meth_duo_grouped)[which(diff_meth_SC_grouped$adjusted_p_UC<=fdr),]
          DUO_UC<-stat_hits[which(abs(stat_hits$db_ctrl_UC)>=dB),]
          nrow(DUO_UC) # 102
          
          # DUO CD
          stat_hits<-as.data.frame(diff_meth_DUO_UCCTRL_organoid)[which(diff_meth_DUO_UCCTRL_organoid$adjusted_p<=fdr),]
          DUO_sig_UCCTRL<-stat_hits[which(abs(stat_hits$delta_beta)>=dB),]
          DUO_sig_UCCTRL$CpG<-rownames(DUO_sig_UCCTRL)
          nrow(DUO_sig_UCCTRL)

#' ### shared hits
intersect(SC_CD$CpG, SC_UC$CpG)
intersect(SC_CD$CpG, TI_sig_UCCTRL$CpG)
intersect(SC_UC$CpG, TI_sig_UCCTRL$CpG)

#' ## NLRC5
TI_sig_UCCTRL[grep("cg07839457",TI_sig_UCCTRL$CpG),]
TI_sig_UCCTRL[grep("cg07862320",TI_sig_UCCTRL$CpG),]

SC_IBD[grep("cg07839457",SC_IBD$CpG),]
SC_IBD[grep("cg07862320",SC_IBD$CpG),]


generic_boxplot(c("cg07839457"))
ggsave(here("DNAm/figs","organoid_CD_SCTIshared_Hit_ThreeBatch.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","organoid_CD_SCTIshared_Hit_ThreeBatch.jpeg"),  width = 7.5, height = 6)

generic_boxplot(c("cg07862320"))






#' ### volcanoes
p.val_max<-min(TI_sig_UCCTRL$p.value[which(TI_sig_UCCTRL$adjusted_p==max(TI_sig_UCCTRL$adjusted_p))])
vol<-makeVolcano(diff_meth_TI_UCCTRL_organoid$p.value, diff_meth_TI_UCCTRL_organoid$delta_beta, 0.05,p.val_max, "Difference Between\nControls + UC and CD", 0.3,11)
ggsave(here("DNAm/figs/jpeg","volcano_TI_CD_organoid_ThreeBatch.jpeg"), vol,  width = 7, height = 6)

p.val_max<-max(SC_CD$p.value_CD[which(SC_CD$adjusted_p_CD==max(SC_CD$adjusted_p_CD))])
vol<-makeVolcano(diff_meth_SC_grouped$p.value_CD, diff_meth_SC_grouped$db_ctrl_CD, 0.05,p.val_max, "Difference Between\nControls and CD", 0.3,11)
ggsave(here("DNAm/figs/jpeg","volcano_SC_CD_organoid_ThreeBatch.jpeg"), vol,  width = 7, height = 6)

p.val_max<-max(SC_UC$p.value_UC[which(SC_UC$adjusted_p_UC==max(SC_UC$adjusted_p_UC))])
vol<-makeVolcano(diff_meth_SC_grouped$p.value_UC, diff_meth_SC_grouped$db_ctrl_UC, 0.05,p.val_max, "Difference Between\nControls and UC", 0.3,11)
ggsave(here("DNAm/figs/jpeg","volcano_SC_UC_organoid_ThreeBatch.jpeg"), vol,  width = 7, height = 6)


#'## R Session Info
sessionInfo()
