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
  beta<-combat_ibd_Beta[cpgs,]
  beta<-melt(beta)
  if(length(cpgs)==1){
    beta$Var2<-rownames(beta)
    beta$Var1<-cpgs}
  plt<-merge(beta, sampleinfo_DNAm, by.x="Var2", by.y="array.id")
  ggplot(plt, aes(diagnosis, value, fill=diagnosis))+
    geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.15, shape=21, size=1,color="black")+
    facet_grid(Var1~sample.site)+fillscale_diagnosis+th+theme_bw()+ylim(0,1)+ylab("DNA Methylation")+xlab("Diagnosis")}






#' ## Load Cell Comp adjusted combatted Pediatric IEC data
load(here("DNAm/data/ibd_adjusted_combatted.Rdata"))


#' ## convert to M values for stats
Mval<-function(beta) log2(beta/(1-beta))
ibd_Mval = apply(combat_ibd_Beta, 1, Mval) # need mvalues for combat
ibd_Mval = as.data.frame(ibd_Mval)
ibd_Mval = t(ibd_Mval)

# only plot a subset fo CpGs for speed
set.seed(1)
mval_melted<- melt(ibd_Mval[sample(nrow(ibd_Mval),10000),])
mval_Plot<-mval_melted[which(!(is.na(mval_melted$value))),]
colnames(mval_Plot)<-c("CpG","ID","Mval")
mval_Plot<-merge(mval_Plot,sampleinfo_DNAm, by.x="ID", by.y="array.id")

ggplot(mval_Plot, aes(Mval, group=as.character(pdata.id), color=as.character(sample.site)))+
  geom_density()+theme_bw()+colscale_sampsite+xlab("DNAm M Value")

ggsave(here("DNAm/figs","Mvalue_distribution_combined.pdf"), w=10, h=5)
ggsave(here("DNAm/figs/jpeg","Mvalue_distribution_combined.jpeg"), w=10, h=5)




#' ## combined meta correlation between potential covariates
meta_categorical <- sampleinfo_DNAm[, c(2,4,5,10,11,13,17)]  # input column numbers in meta that contain categorical variables "case.no"      "diagnosis"    "sample.site"  "inflammation" "sex"          "array.type"   "sample_ID" 
meta_continuous <- as.data.frame(sampleinfo_DNAm[, c(12)] ) # input column numbers in meta that contain continuous variables
colnames(meta_continuous) <- c("Age")

Meta_correlation<-lapply(1:ncol(meta_continuous), function(x) as.numeric(meta_continuous[,x]))
Meta_correlation<-as.data.frame(do.call(cbind, Meta_correlation))
colnames(Meta_correlation)<-colnames(meta_continuous)

aov_meta <- sapply(1:ncol(meta_categorical), function(cat) {
  summary(aov(meta_continuous$Age ~ meta_categorical[, cat]))[[1]]$"Pr(>F)"[1]}
)

suppressWarnings(chi_meta<-lapply(1:ncol(meta_categorical), function(cat1) {
  sapply(1:ncol(meta_categorical), function(cat2){
    chisq.test(table(meta_categorical[,cat1], meta_categorical[,cat2]))$p.value})
}))

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
  fillscale_pval

ggsave(here("DNAm/figs","primary_meta_correlation_both_arrays.pdf"), width = 6, height = 5)
ggsave(here("DNAm/figs/jpeg","meta_correlation_both_arrays.jpeg"), width = 6, height = 5)

#' ### proportion of diagnosis and sample site on each array
pro_array<-as.data.frame.matrix(table(sampleinfo_DNAm$diagnosis, sampleinfo_DNAm$array.type))
pro_array$EPIC<-(pro_array$EPIC/sum(pro_array$EPIC))*100
pro_array$K450<-(pro_array$K450/sum(pro_array$K450))*100
pro_array$diagnosis<-rownames(pro_array)
pro_array<-melt(pro_array)

ggplot(pro_array, aes(diagnosis,value, fill=variable))+geom_bar(stat="identity", position="dodge")+theme_bw()

pro_arraysite<-as.data.frame.matrix(table(sampleinfo_DNAm$sample.site, sampleinfo_DNAm$array.type))
pro_arraysite$EPIC<-(pro_arraysite$EPIC/sum(pro_arraysite$EPIC))*100
pro_arraysite$K450<-(pro_arraysite$K450/sum(pro_arraysite$K450))*100
pro_arraysite$diagnosis<-rownames(pro_arraysite)
pro_arraysite<-melt(pro_arraysite)

ggplot(pro_arraysite, aes(diagnosis,value, fill=variable))+geom_bar(stat="identity", position="dodge")+theme_bw()

#' ### Inflammation and sample site
pro_inflamed<-as.data.frame.matrix(table(sampleinfo_DNAm$sample.site, sampleinfo_DNAm$inflammation))
pro_inflamed$total<-rowSums(pro_inflamed)
pro_inflamed$perc_inflamed<-(pro_inflamed$Y/pro_inflamed$total)*100
pro_inflamed

#' ### Inflammation and diagnosis
pro_inflamed_dia<-as.data.frame.matrix(table(sampleinfo_DNAm$diagnosis, sampleinfo_DNAm$inflammation))
pro_inflamed_dia$total<-rowSums(pro_inflamed_dia)
pro_inflamed_dia$perc_inflamed<-(pro_inflamed_dia$Y/pro_inflamed_dia$total)*100
pro_inflamed_dia





#' ## remove rescopes 
#' They will be analysed in a separate analysis
sampleinfo_og<-sampleinfo_DNAm[which(is.na(sampleinfo_DNAm$passage.or.rescope.no)),]
ibd_Mval_og<-ibd_Mval[,which(colnames(ibd_Mval)%in%sampleinfo_og$array.id)]
ibd_combo_og<-combat_ibd_Beta[,which(colnames(combat_ibd_Beta)%in%sampleinfo_og$array.id)]
identical(colnames(ibd_Mval_og),sampleinfo_og$array.id)

#'## IBD vs controls in each tissue
sampleinfo_og$diagnosis_simple<-as.factor(sampleinfo_og$diagnosis)
levels(sampleinfo_og$diagnosis_simple)<-c("IBD","Control","IBD","IBD","IBD")



#' ## subset  to individual sampling sites
sampleinfo_og_TI<-sampleinfo_og[which(sampleinfo_og$sample.site=="TI"),]
ibd_Mval_og_TI<-ibd_Mval[,which(colnames(ibd_Mval)%in%sampleinfo_og_TI$array.id)]
ibd_combo_og_TI<-combat_ibd_Beta[,which(colnames(combat_ibd_Beta)%in%sampleinfo_og_TI$array.id)]
identical(colnames(ibd_Mval_og_TI),sampleinfo_og_TI$array.id)

sampleinfo_og_AC<-sampleinfo_og[which(sampleinfo_og$sample.site=="AC"),]
ibd_Mval_og_AC<-ibd_Mval[,which(colnames(ibd_Mval)%in%sampleinfo_og_AC$array.id)]
ibd_combo_og_AC<-combat_ibd_Beta[,which(colnames(combat_ibd_Beta)%in%sampleinfo_og_AC$array.id)]
identical(colnames(ibd_Mval_og_AC),sampleinfo_og_AC$array.id)

sampleinfo_og_SC<-sampleinfo_og[which(sampleinfo_og$sample.site=="SC"),]
ibd_Mval_og_SC<-ibd_Mval[,which(colnames(ibd_Mval)%in%sampleinfo_og_SC$array.id)]
ibd_combo_og_SC<-combat_ibd_Beta[,which(colnames(combat_ibd_Beta)%in%sampleinfo_og_SC$array.id)]
identical(colnames(ibd_Mval_og_SC),sampleinfo_og_SC$array.id)



#######################
#'##  IBD vs controls 
#' only relevant in SC
#######################
DMA<-function(mval,beta,meta){
  print(paste("Testing in ", ncol(beta), " individuals and only ", unique(meta$sample.site), " samples", sep=""))
  # p value
  mod<-model.matrix(~as.factor(diagnosis_simple)+as.factor(inflammation)+as.factor(sex)+age, data=meta)
  fit <- lmFit(mval, mod)
  ebfit <- eBayes(fit)
  
  # covariate adjusted beta values
  avebeta.lm<-apply(beta, 1, function(x){
    lm(x~as.factor(inflammation)+as.factor(sex)+age, data=meta)})
  residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
  colnames(residuals)<-colnames(beta)
  adj.residuals<-residuals+matrix(apply(beta, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
  
  delta_beta<-sapply(1:nrow(adj.residuals), function(x) {
    group_means<-tapply(adj.residuals[x,], meta$diagnosis_simple, mean)
    group_means[1]-group_means[2]
  })
  
  limmares <- data.frame(p.value=ebfit$p.value[,"as.factor(diagnosis_simple)Control"], delta_beta=delta_beta)
  limmares}


diff_meth_TI<-DMA(ibd_Mval_og_TI,ibd_combo_og_TI, sampleinfo_og_TI)
diff_meth_AC<-DMA(ibd_Mval_og_AC,ibd_combo_og_AC, sampleinfo_og_AC)
diff_meth_SC<-DMA(ibd_Mval_og_SC,ibd_combo_og_SC, sampleinfo_og_SC)

diff_meth_TI$adjusted_p<-p.adjust(diff_meth_TI$p.value, method = "fdr", n = nrow(diff_meth_TI))
diff_meth_AC$adjusted_p<-p.adjust(diff_meth_AC$p.value, method = "fdr", n = nrow(diff_meth_AC))
diff_meth_SC$adjusted_p<-p.adjust(diff_meth_SC$p.value, method = "fdr", n = nrow(diff_meth_SC))

save(diff_meth_TI, diff_meth_AC,diff_meth_SC, file=here("DNAm/data","diff_DNAm_simple_casectrl.RData"))




#######################
#' ## Differential DNAm in disease subset specific grouping
#######################
DMA_grouped<-function(mval,beta,meta){
  print(paste("Testing in ", ncol(beta), " individuals and only ", unique(meta$sample.site), " samples", sep=""))
  # p value
  meta$diagnosis_grouped<-as.factor(meta$diagnosis_grouped)
  meta$inflammation<-as.factor(meta$inflammation)
  meta$sex<-as.factor(meta$sex)
  
  mod<-model.matrix(~0+diagnosis_grouped+inflammation+sex+age, data=meta)
  fit <- lmFit(mval, mod)
  #ebfit <- eBayes(fit)
  
  
 # construct the contrast matrix
  contrastMatrix <- makeContrasts(
    ctrl_CD = diagnosis_groupedControl - diagnosis_groupedCD,
    ctrl_UC = diagnosis_groupedControl - diagnosis_groupedUC,
    CD_UC = diagnosis_groupedCD - diagnosis_groupedUC,
    levels = mod
  )
  
  contrastFit<-contrasts.fit(fit, contrastMatrix)
  contrastFitEb <- eBayes(contrastFit) 
  contrastCpG_CD <- topTable(contrastFitEb, coef="ctrl_CD",number=Inf)
  contrastCpG_UC <- topTable(contrastFitEb, coef="ctrl_UC",number=Inf)
  contrastCpG_CDUC <- topTable(contrastFitEb, coef="CD_UC",number=Inf)
  
  # covariate adjusted beta values
  avebeta.lm<-apply(beta, 1, function(x){
    lm(x~as.factor(inflammation)+as.factor(sex)+age, data=meta)})
  residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
  colnames(residuals)<-colnames(beta)
  adj.residuals<-residuals+matrix(apply(beta, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
  
  delta_beta<-lapply(1:nrow(adj.residuals), function(x) {
    group_means<-tapply(adj.residuals[x,], meta$diagnosis_grouped, mean)
    c(group_means[1]-group_means[2],group_means[1]-group_means[3], group_means[3]-group_means[2])
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

diff_meth_TI_grouped<-DMA_grouped(ibd_Mval_og_TI,ibd_combo_og_TI, sampleinfo_og_TI)
diff_meth_AC_grouped<-DMA_grouped(ibd_Mval_og_AC,ibd_combo_og_AC, sampleinfo_og_AC)
diff_meth_SC_grouped<-DMA_grouped(ibd_Mval_og_SC,ibd_combo_og_SC, sampleinfo_og_SC)

save(diff_meth_TI_grouped, diff_meth_AC_grouped,diff_meth_SC_grouped, file=here("DNAm/data/diff_DNAm_grouped_casectrl.RData"))



####################
#' # Grouping UC with controls in TI makes biological sense
####################

# IBD vs Ctrl (done above)
# UC vs Ctrl in SC (done above)
# CD vs Ctrl in SC (done above)
# CD vs UC+Ctrl in TI (done below)

mval<-ibd_Mval_og_TI
beta<-ibd_combo_og_TI
meta<-sampleinfo_og_TI

print(paste("Testing in ", ncol(beta), " individuals and only ", unique(meta$sample.site), " samples", sep=""))
# p value
meta$diagnosis_grouped<-as.factor(meta$diagnosis_grouped)
meta$diagnosis_TI_split<-meta$diagnosis_grouped
levels(meta$diagnosis_TI_split)<-c("Control_UC","CD","Control_UC")

meta$inflammation<-as.factor(meta$inflammation)
meta$sex<-as.factor(meta$sex)

mod<-model.matrix(~diagnosis_TI_split+inflammation+sex+age, data=meta)
fit <- lmFit(mval, mod)
ebfit <- eBayes(fit)

avebeta.lm<-apply(beta, 1, function(x){
  lm(x~as.factor(inflammation)+as.factor(sex)+age, data=meta)})
residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
colnames(residuals)<-colnames(beta)
adj.residuals<-residuals+matrix(apply(beta, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))

delta_beta<-sapply(1:nrow(adj.residuals), function(x) {
  group_means<-tapply(adj.residuals[x,], meta$diagnosis_TI_split, mean)
  group_means[1]-group_means[2] # Control_UC - CD
})

diff_meth_TI_UCCTRL <- data.frame(p.value=ebfit$p.value[,"diagnosis_TI_splitCD"], delta_beta=delta_beta)
save(diff_meth_TI_UCCTRL, file=here("DNAm/data/diff_DNAm_TI_UCCTRL.RData"))




#' ### load all four models
# load(file=here("DNAm/data/diff_DNAm_simple_casectrl.RData"))
# load(file=here("DNAm/data/diff_DNAm_grouped_casectrl.RData"))
# load(file=here("DNAm/data/diff_DNAm_TI_UCCTRL.RData"))

diff_meth_TI_UCCTRL$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL))


#' ## Significant CpGs numbers
dB=0.05
fdr=0.05


# SC IBD
stat_hits<-as.data.frame(diff_meth_SC)[which(diff_meth_SC$adjusted_p<=fdr),]
SC_IBD<-stat_hits[which(abs(stat_hits$delta_beta)>=dB),]
nrow(SC_IBD)
SC_IBD$CpG<-rownames(SC_IBD)


# SC CD
stat_hits<-as.data.frame(diff_meth_SC_grouped)[which(diff_meth_SC_grouped$adjusted_p_CD<=fdr),]
SC_CD<-stat_hits[which(abs(stat_hits$db_ctrl_CD)>=dB),]
nrow(SC_CD)

# SC UC
stat_hits<-as.data.frame(diff_meth_SC_grouped)[which(diff_meth_SC_grouped$adjusted_p_UC<=fdr),]
SC_UC<-stat_hits[which(abs(stat_hits$db_ctrl_UC)>=dB),]
nrow(SC_UC)

# TI CD
stat_hits<-as.data.frame(diff_meth_TI_UCCTRL)[which(diff_meth_TI_UCCTRL$adjusted_p<=fdr),]
TI_sig_UCCTRL<-stat_hits[which(abs(stat_hits$delta_beta)>=dB),]
TI_sig_UCCTRL$CpG<-rownames(TI_sig_UCCTRL)
nrow(TI_sig_UCCTRL)

#' ### no shared hits
intersect(SC_CD$CpG, SC_UC$CpG)
intersect(SC_CD$CpG, TI_sig_UCCTRL$CpG)
intersect(SC_UC$CpG, TI_sig_UCCTRL$CpG)

#' ## NLRC5
TI_sig_UCCTRL[grep("cg07839457",TI_sig_UCCTRL$CpG),]
#' not tested on 450K so filtered as not available in whole cohort
TI_sig_UCCTRL[grep("cg07862320",TI_sig_UCCTRL$CpG),]




generic_boxplot(c("cg02004044","cg01098981"))
ggsave(here("DNAm/figs","EWAS_random_Hit.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","EWAS_random_Hit.jpeg"), width = 7.5, height = 6)

#'  TI UCCTRL
generic_boxplot(rownames(TI_sig_UCCTRL)[rev(order(abs(TI_sig_UCCTRL$delta_beta)))][1:4])



#'### volcanoes
#'sub in the unadjusted p value max at FDR 0.05, will be slightly different for each comparison
p.val_max<-min(TI_sig_UCCTRL$p.value[which(TI_sig_UCCTRL$adjusted_p==max(TI_sig_UCCTRL$adjusted_p))])
vol<-makeVolcano(diff_meth_TI_UCCTRL$p.value, diff_meth_TI_UCCTRL$delta_beta, 0.05,p.val_max, "Difference Between\nControls+UC and CD", 0.3, 8.5)
ggsave(here("DNAm/figs/jpeg","volcano_TI_CD.jpeg"),vol,  width = 8, height = 6)

p.val_max<-max(SC_CD$p.value_CD[which(SC_CD$adjusted_p_CD==max(SC_CD$adjusted_p_CD))])
vol<-makeVolcano(diff_meth_SC_grouped$p.value_CD, diff_meth_SC_grouped$db_ctrl_CD, 0.05,p.val_max, "Difference Between\nControls and CD", 0.3, 8.5)
ggsave(here("DNAm/figs/jpeg","volcano_SC_CD.jpeg"),vol,  width = 8, height = 6)

p.val_max<-max(SC_UC$p.value_UC[which(SC_UC$adjusted_p_UC==max(SC_UC$adjusted_p_UC))])
vol<-makeVolcano(diff_meth_SC_grouped$p.value_UC, diff_meth_SC_grouped$db_ctrl_UC, 0.05,9e-07, "Difference Between\nControls and UC", 0.3, 8.5)
ggsave(here("DNAm/figs/jpeg","volcano_SC_UC.jpeg"),vol,  width = 8, height = 6)



#' ## Compare to howell differential CpGs

#' ### CD

TI_Ctrl_CD<-read.table(here("DNAm/data/public/howell_data/Howell_Gastroenterology_dffDNAm_TI_CTRLCD.csv"), skip=1, sep=",", header=T)
length(intersect(TI_Ctrl_CD$X, TI_sig_UCCTRL$CpG))


TI_sig_UCCTRL_howell<-TI_sig_UCCTRL[which(TI_sig_UCCTRL$CpG%in%TI_Ctrl_CD$X),]
TI_Ctrl_CD<-TI_Ctrl_CD[which(TI_Ctrl_CD$X%in%TI_sig_UCCTRL_howell$CpG),]
TI_Ctrl_CD<-TI_Ctrl_CD[match(TI_sig_UCCTRL_howell$CpG,TI_Ctrl_CD$X),]
identical(TI_Ctrl_CD$X, TI_sig_UCCTRL_howell$CpG)

cor(TI_Ctrl_CD$P.Value, TI_sig_UCCTRL_howell$p.value)
cor(TI_Ctrl_CD$logFC, TI_sig_UCCTRL_howell$delta_beta)



SC_Ctrl_CD<-read.table(here("DNAm/data/public/howell_data/Howell_Gastroenterology_dffDNAm_SC_CTRLCD.csv"), skip=1, sep=",", header=T)
length(intersect(SC_Ctrl_CD$X, SC_CD$CpG))

SC_sig_CD_howell<-SC_CD[which(SC_CD$CpG%in%SC_Ctrl_CD$X),]
SC_Ctrl_CD<-SC_Ctrl_CD[which(SC_Ctrl_CD$X%in%SC_sig_CD_howell$CpG),]
SC_Ctrl_CD<-SC_Ctrl_CD[match(SC_sig_CD_howell$CpG,SC_Ctrl_CD$X),]
identical(SC_Ctrl_CD$X, SC_sig_CD_howell$CpG)

cor(SC_Ctrl_CD$P.Value, SC_sig_CD_howell$p.value_CD)
cor(SC_Ctrl_CD$logFC, SC_sig_CD_howell$db_ctrl_CD)



#' ### UC

SC_Ctrl_UC<-read.table(here("DNAm/data/public/howell_data/Howell_Gastroenterology_dffDNAm_SC_CTRLUC.csv"), skip=1, sep=",", header=T)
length(intersect(SC_Ctrl_UC$X, SC_CD$CpG))

SC_sig_UC_howell<-SC_UC[which(SC_UC$CpG%in%SC_Ctrl_UC$X),]
SC_Ctrl_UC<-SC_Ctrl_UC[which(SC_Ctrl_UC$X%in%SC_sig_UC_howell$CpG),]
SC_Ctrl_UC<-SC_Ctrl_UC[match(SC_sig_UC_howell$CpG,SC_Ctrl_UC$X),]
identical(SC_Ctrl_UC$X, SC_sig_UC_howell$CpG)

print(c(SC_Ctrl_UC$P.Value, SC_sig_UC_howell$p.value_UC))
print(c(SC_Ctrl_UC$logFC, SC_sig_UC_howell$db_ctrl_UC))




#'## R Session Info
sessionInfo()



