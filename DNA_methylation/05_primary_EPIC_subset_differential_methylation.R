#'---
#'title: Primary differenital DNAm in CpGs only measured on EPIC
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ## Cell adjust CpGs not done in the main dataset
#' Those were only CpGs on both EPIC and 450K 

#' Adapted from the estimateCellCounts function from minfi and ECC7 script by Dr. Meaghan Jones

#' ### Load Libraries
#' require(quadprog) only needed fro more than 2 cell types
suppressMessages(library(reshape))
library(genefilter)
suppressMessages(library(matrixStats))
library(limma)
library(ggplot2)
library(gridExtra)
library(quadprog)
suppressMessages(library(sva))
suppressMessages(library(pamr))
library(limma)
library(here)

options(stringsAsFactors = FALSE)

source(here("DNAm/scripts/","04_minfi_estcellcounts_functions.R"))
source(here("general_functions","00_pretty_plots.R"))
source(here("general_functions","00_Heat_scree_plot_generic.R"))



#' Load in EPIC data with is functionally normalized but not cell comp adjusted
load(here("../data/ibd_beta_botharrays_funnorm_filtered.RData"))
dim(ibd_beta_epic)

#' Load The cell comp counts for all samples to adjust by the same amount at these epic only CpGs as well
load(here("DNAm/data/","DNAm_organoid_counts.RData"))


sampleinfo_DNAm<-sampleinfo_DNAm_countsConstrained
sampleinfo_DNAm$diagnosis_grouped<-as.factor(sampleinfo_DNAm$diagnosis)
levels(sampleinfo_DNAm$diagnosis_grouped)<-c("CD","Control","CD","UC","UC")
sampleinfo_DNAm$diagnosis_grouped<-factor(sampleinfo_DNAm$diagnosis_grouped, levels=c("Control", "CD", "UC"))
sampleinfo_DNAm$sample.site_grouped<-as.factor(sampleinfo_DNAm$sample.site)
levels(sampleinfo_DNAm$sample.site_grouped)<-c("colon","colon","ileum")


#'## Just EPIC samples
sampleinfo_DNAm<-sampleinfo_DNAm[which(sampleinfo_DNAm$array.type=="EPIC"),]
sampleinfo_DNAm<-sampleinfo_DNAm[match(colnames(ibd_beta_epic),sampleinfo_DNAm$array.id),]
identical(colnames(ibd_beta_epic), sampleinfo_DNAm$array.id)



### Cell composition adjustment

avebeta.lm<-apply(ibd_beta_epic, 1, function(x){
  lm(x~sampleinfo_DNAm$organoid)
})
residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
colnames(residuals)<-colnames(ibd_beta_epic)
ibd_epic_adjusted<-residuals+matrix(apply(ibd_beta_epic, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))

save(ibd_epic_adjusted, file=paste(here("DNAm/data/"),"ibd_adjusted_betas_EPIC_only.RData", sep=""))




#' # Differential methylation
#' ### For efficiency only test CpGs not tested on both arrays
load(file=here("DNAm/data/diff_DNAm_grouped_casectrl.RData"))
ibd_beta_epic<-ibd_epic_adjusted[which(!(rownames(ibd_epic_adjusted)%in%diff_meth_TI_grouped$CpG)),]
dim(ibd_beta_epic)

# convert to M values
Mval<-function(beta) log2(beta/(1-beta))
ibd_Mval = suppressWarnings(apply(ibd_beta_epic, 1, Mval)) # need mvalues for combat
ibd_Mval = as.data.frame(ibd_Mval)
ibd_Mval = t(ibd_Mval)


# remove rescopes 
sampleinfo_og<-sampleinfo_DNAm[which(is.na(sampleinfo_DNAm$passage.or.rescope.no)),]
ibd_Mval_og<-ibd_Mval[,which(colnames(ibd_Mval)%in%sampleinfo_og$array.id)]
ibd_combo_og<-ibd_beta_epic[,which(colnames(ibd_beta_epic)%in%sampleinfo_og$array.id)]
identical(colnames(ibd_Mval_og),sampleinfo_og$array.id)

#'## IBD vs controls in each tissue
sampleinfo_og$diagnosis_simple<-as.factor(sampleinfo_og$diagnosis)
levels(sampleinfo_og$diagnosis_simple)<-c("IBD","Control","IBD","IBD","IBD")


#subset tissue
sampleinfo_og_TI<-sampleinfo_og[which(sampleinfo_og$sample.site=="TI"),]
ibd_Mval_og_TI<-ibd_Mval[,which(colnames(ibd_Mval)%in%sampleinfo_og_TI$array.id)]
ibd_combo_og_TI<-ibd_combo_og[,which(colnames(ibd_combo_og)%in%sampleinfo_og_TI$array.id)]
identical(colnames(ibd_Mval_og_TI),sampleinfo_og_TI$array.id)

sampleinfo_og_SC<-sampleinfo_og[which(sampleinfo_og$sample.site=="SC"),]
ibd_Mval_og_SC<-ibd_Mval[,which(colnames(ibd_Mval)%in%sampleinfo_og_SC$array.id)]
ibd_combo_og_SC<-ibd_combo_og[,which(colnames(ibd_combo_og)%in%sampleinfo_og_SC$array.id)]
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
diff_meth_SC<-DMA(ibd_Mval_og_SC,ibd_combo_og_SC, sampleinfo_og_SC)

diff_meth_TI$adjusted_p<-p.adjust(diff_meth_TI$p.value, method = "fdr", n = nrow(diff_meth_TI))
diff_meth_SC$adjusted_p<-p.adjust(diff_meth_SC$p.value, method = "fdr", n = nrow(diff_meth_SC))






#######################
#' ### disease subset specific grouping
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


diff_meth_SC_grouped<-DMA_grouped(ibd_Mval_og_SC,ibd_combo_og_SC, sampleinfo_og_SC)


####################
## More biologically driven model for TI Matthais suggested
####################

# UC vs Ctrl in SC (done already)
# CD vs Ctrl in SC (done already)
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
save(diff_meth_TI_UCCTRL,diff_meth_SC_grouped,diff_meth_TI,diff_meth_SC,  file=paste(here("DNAm/data/"),"diff_DNAm_EPIC_only.RData", sep=""))



#' ## Look at CpGs
diff_meth_TI_UCCTRL$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL))


## Significant CpGs exploration
stat_hits<-as.data.frame(diff_meth_SC_grouped)[which(diff_meth_SC_grouped$adjusted_p_CD<=0.05),]
SC_CD<-stat_hits[which(abs(stat_hits$db_ctrl_CD)>=0.05),]
nrow(SC_CD)
stat_hits<-as.data.frame(diff_meth_SC_grouped)[which(diff_meth_SC_grouped$adjusted_p_UC<=0.05),]
SC_UC<-stat_hits[which(abs(stat_hits$db_ctrl_UC)>=0.05),]
nrow(SC_UC)

stat_hits<-as.data.frame(diff_meth_TI_UCCTRL)[which(diff_meth_TI_UCCTRL$adjusted_p<=0.05),]
TI_sig_UCCTRL<-stat_hits[which(abs(stat_hits$delta_beta)>=0.05),]
TI_sig_UCCTRL$CpG<-rownames(TI_sig_UCCTRL)
nrow(TI_sig_UCCTRL)

intersect(SC_CD$CpG, SC_UC$CpG)
intersect(SC_CD$CpG, TI_sig_UCCTRL$CpG)
intersect(SC_UC$CpG, TI_sig_UCCTRL$CpG)


#'### from gene_level analysis
#' TAP1
diff_meth_TI_UCCTRL[grep("cg11706729",rownames(diff_meth_TI_UCCTRL)),]
diff_meth_TI_UCCTRL[grep("cg02756056",rownames(diff_meth_TI_UCCTRL)),]

#' NLRC5
diff_meth_TI_UCCTRL[grep("cg07839457",rownames(diff_meth_TI_UCCTRL)),]
diff_meth_TI_UCCTRL[grep("cg07862320",rownames(diff_meth_TI_UCCTRL)),]

#' B2M
diff_meth_TI_UCCTRL[grep("cg27537252",rownames(diff_meth_TI_UCCTRL)),]
diff_meth_TI_UCCTRL[grep("cg02988397",rownames(diff_meth_TI_UCCTRL)),]


