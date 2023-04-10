#'---
#'title: Differential DNAm Rescopes
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(limma)
library(Metrics)
library(here)

options(stringsAsFactors = FALSE)

#' functions
source(here("general_functions","00_pretty_plots.R"))
source(here("general_functions","00_Heat_scree_plot_generic.R"))


load(here("DNAm/data","ibd_adjusted_combatted.Rdata"))
# load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/ibd_adjusted_combatted.Rdata"))


# convert to M values
Mval<-function(beta) log2(beta/(1-beta))
ibd_Mval = apply(combat_ibd_Beta, 1, Mval) # need mvalues for combat
ibd_Mval = as.data.frame(ibd_Mval)
ibd_Mval = t(ibd_Mval)




#'## rescope RMSE
sampleinfo_rescope<-sampleinfo_DNAm[which(sampleinfo_DNAm$passage.or.rescope.no=="RE1"),]
rescope_samples<-unique(sampleinfo_rescope$sample_ID)
sampleinfo_rescope<-sampleinfo_DNAm[which(sampleinfo_DNAm$sample_ID%in%rescope_samples),]

rescope_rmse<-lapply(1:length(rescope_samples), function(x){
  samples<-sampleinfo_rescope$array.id[which(sampleinfo_rescope$sample_ID==rescope_samples[x])]
  if(length(samples)==1){}else{
    df_rmse<-combat_ibd_Beta[,samples]
    data.frame(RMSE=rmse(df_rmse[,1], df_rmse[,2]), correlation=cor(df_rmse[,1], df_rmse[,2], method="spearman"),sample=rescope_samples[x])}
})
rescope_rmse<-do.call(rbind, rescope_rmse)

ggplot(rescope_rmse, aes(RMSE))+geom_density()
ggsave(here("DNAm/figs","rescope_RMSE.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","rescope_RMSE.jpeg"),  width = 7.5, height = 6)




#' # Validation of differential methylation
# take only rescopes n=47
sampleinfo_re<-sampleinfo_DNAm[which(!(is.na(sampleinfo_DNAm$passage.or.rescope.no))),]
ibd_Mval_re<-ibd_Mval[,which(colnames(ibd_Mval)%in%sampleinfo_re$array.id)]
ibd_combo_re<-combat_ibd_Beta[,which(colnames(combat_ibd_Beta)%in%sampleinfo_re$array.id)]
identical(colnames(ibd_combo_re),sampleinfo_re$array.id)

tapply(sampleinfo_re$array.id, list(sampleinfo_re$sample.site, sampleinfo_re$diagnosis_grouped), length)

# No controls in rescope samples, so can do differential analysis with a subset of controls or just show DNAm consistent in rescopes at DIff CpGs
load(file=here("DNAm/data","diff_DNAm_grouped_casectrl.RData"))
#load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/diff_DNAm_grouped_casectrl.RData"))

load(file=here("DNAm/data","diff_DNAm_TI_UCCTRL.RData"))
#load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/diff_DNAm_TI_UCCTRL.RData"))
diff_meth_TI_UCCTRL$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL))


## Significant CpGs exploration
# CD
stat_hits<-as.data.frame(diff_meth_SC_grouped)[which(diff_meth_SC_grouped$adjusted_p_CD<=0.05),]
SC_CD<-stat_hits[which(abs(stat_hits$db_ctrl_CD)>=0.05),]

# UC
stat_hits<-as.data.frame(diff_meth_TI_UCCTRL)[which(diff_meth_TI_UCCTRL$adjusted_p<=0.05),]
TI_sig_UCCTRL<-stat_hits[which(abs(stat_hits$delta_beta)>=0.05),]
TI_sig_UCCTRL$CpG<-rownames(TI_sig_UCCTRL)

stat_hits<-as.data.frame(diff_meth_SC_grouped)[which(diff_meth_SC_grouped$adjusted_p_UC<=0.05),]
SC_UC<-stat_hits[which(abs(stat_hits$db_ctrl_UC)>=0.05),]



#' ## bland altman and correlation at differential CpGs
sampleinfo_DNAm$rescope_plt<-sampleinfo_DNAm$passage.or.rescope.no
sampleinfo_DNAm$rescope_plt[is.na(sampleinfo_DNAm$rescope_plt)]<-"OG"
sampleinfo_DNAm$rescope_plt<-as.factor(sampleinfo_DNAm$rescope_plt)
levels(sampleinfo_DNAm$rescope_plt)<-c("Original","Rescope")



plt_of_CpG_val<-function(CpG, tissue){
  
  smp_info_plt<-sampleinfo_DNAm[which(sampleinfo_DNAm$sample.site==tissue),]
  smp_info_plt$betas<-combat_ibd_Beta[CpG,which(sampleinfo_DNAm$sample.site==tissue)]
  smp_info_plt$diagnosis_grouped<-factor(smp_info_plt$diagnosis_grouped, levels=c("Control","UC","CD"))
  
  
  box<-ggplot(smp_info_plt, aes(rescope_plt, betas, fill=interaction(diagnosis_grouped, rescope_plt)))+
    geom_boxplot(outlier.shape=NA)+geom_line(aes(group=sample_ID), color="lightgrey")+
    geom_point(position = position_dodge(width=0.75), shape=21, size=2)+th+theme_bw()+facet_wrap(~diagnosis_grouped)+ylim(0,1)+
    scale_fill_manual(values=c("lightgrey", "darkgoldenrod1", "dodgerblue3", "#ffe39f","#b9d5f0"), guide=F)+
    xlab("Sample Time Point")+ylab("DNAm Beta")+ggtitle(CpG)
  
  sampleinfo_re_ID<-sampleinfo_DNAm[which(!(is.na(sampleinfo_DNAm$passage.or.rescope.no))),"sample_ID"]
  sampleinfo_DNAm_reog<-sampleinfo_DNAm[which(sampleinfo_DNAm$sample_ID%in%sampleinfo_re_ID),]
  
  ibd_combo_reog<-combat_ibd_Beta[,which(colnames(combat_ibd_Beta)%in%sampleinfo_DNAm_reog$array.id)]
  identical(colnames(ibd_combo_reog),sampleinfo_DNAm_reog$array.id)
  
  smp_info_plt<-sampleinfo_DNAm_reog[which(sampleinfo_DNAm_reog$sample.site==tissue),]
  smp_info_plt$betas<-ibd_combo_reog[CpG,which(sampleinfo_DNAm_reog$sample.site==tissue)]
  
  og<-smp_info_plt[which(is.na(smp_info_plt$passage.or.rescope.no)),]
  re<-smp_info_plt[which(!(is.na(smp_info_plt$passage.or.rescope.no))),]
  re<-re[which(re$case.no%in%og$case.no),]
  re<-re[order(re$case.no),]
  og<-og[order(og$case.no),]
  identical(re$case.no, og$case.no)
  
  plt<-data.frame(case.no=re$case.no, re_beta=re$betas, og_beta=og$betas)
  plt$mean<-rowMeans(plt[,2:3])
  plt$diff<-plt$og_beta-plt$re_beta
  
  corcpg=round(cor(re$betas,og$betas), 2)
  rmsecpg=round(rmse(re$betas,og$betas), 2)
  
  plt<-merge(plt, smp_info_plt[,c(2,21)], by="case.no")
  
  cor<-ggplot(plt, aes(og_beta, re_beta))+geom_point(aes(fill=diagnosis_grouped),shape=21, color="black", size=2)+fillscale_diagnosis+
    stat_smooth(method="lm", se=F, color="black")+ylim(0,1)+xlim(0,1)+theme_bw()+annotate("text", x=0.15, y=0.8, label=paste("R = ",corcpg, sep=""))+
    theme(legend.position = "none")+xlab("Original Sample DNAm")+ylab("Rescope Sample DNAm")
  bland<-ggplot(plt, aes(mean, diff))+geom_point(aes(fill=diagnosis_grouped),shape=21,  color="black",size=2)+fillscale_diagnosis+
    geom_hline(yintercept=0)+theme_bw()+annotate("text", x=0.45, y=0.25, label=paste("RMSE = ", rmsecpg, sep=""))+
    theme(legend.position = "none")+xlab("Mean DNAm of\nOriginal and Rescope")+ylab("DNAm Difference Between\nOriginal and Rescope")
  grid.arrange(box, arrangeGrob(cor, bland, ncol=2), heights=c(1,0.75))
}


plt_of_CpG_val("cg06708937", "SC")
plt_of_CpG_val("cg22072829", "SC")

#' #' ### TAP1 
#' #' CpG cg24898914 cg25384897 in organoid and primary
#' cpg1<-plt_of_CpG_val("cg24898914", "TI")
#' cpg2<-plt_of_CpG_val("cg25384897", "TI")
#' grid.arrange(cpg1,cpg2, ncol=2)
#' ggsave(file=here("DNAm/figs","TAP1_rescopes.pdf"),grid.arrange(cpg1,cpg2, ncol=2), width=12, height=6)
#' ggsave(file=here("DNAm/figs/jpeg","TAP1_rescopes.jpeg"),grid.arrange(cpg1,cpg2, ncol=2), width=12, height=6)


#' ### TAP1 
#' CpG cg02756056 in primary (cg11706729 was sig but only in organoids)
cpg1<-plt_of_CpG_val("cg02756056", "TI")
cpg1
ggsave(file=here("DNAm/figs","TAP1_rescopes.pdf"),cpg1, width=6, height=6)
ggsave(file=here("DNAm/figs/jpeg","TAP1_rescopes.jpeg"),cpg1, width=6, height=6)



#' ### NLRC5 
#' CpG cg07839457 (cg07862320 only on EPIC) in organoid and primary
cpg1<-plt_of_CpG_val("cg07839457", "TI")
cpg1
ggsave(file=here("DNAm/figs","NLRC5_rescopes.pdf"),cpg1, width=6, height=6)
ggsave(file=here("DNAm/figs/jpeg","NLRC5_rescopes.jpeg"),cpg1, width=6, height=6)


#' ### MUC1 
#' CpG cg24512973 cg18804777" in organoid and primary 
cpg1<-plt_of_CpG_val("cg24512973", "TI")
cpg2<-plt_of_CpG_val("cg18804777", "TI")
grid.arrange(cpg1,cpg2, ncol=2)
ggsave(file=here("DNAm/figs","MUC1_rescopes.pdf"),grid.arrange(cpg1,cpg2, ncol=2), width=12, height=6)
ggsave(file=here("DNAm/figs/jpeg","MUC1_rescopes.jpeg"),grid.arrange(cpg1,cpg2, ncol=2), width=12, height=6)



### rescope for presentation

plt_of_CpG_val_less<-function(CpG, tissue){
  smp_info_plt<-sampleinfo_DNAm[which(sampleinfo_DNAm$sample.site==tissue),]
  smp_info_plt$betas<-combat_ibd_Beta[CpG,which(sampleinfo_DNAm$sample.site==tissue)]
  smp_info_plt$diagnosis_grouped<-factor(smp_info_plt$diagnosis_grouped, levels=c("Control","UC","CD"))
  box<-ggplot(smp_info_plt, aes(rescope_plt, betas, fill=interaction(diagnosis_grouped, rescope_plt)))+
    geom_boxplot(outlier.shape=NA)+
    geom_point(position = position_dodge(width=0.75), shape=21, size=2)+th_present+theme_bw()+facet_wrap(~diagnosis_grouped)+ylim(0,1)+
    scale_fill_manual(values=c("lightgrey", "darkgoldenrod1", "dodgerblue3", "#ffe39f","#b9d5f0"), guide=F)+
    xlab("Sample Time Point")+ylab("DNAm Beta")+ggtitle(CpG)+
    theme(legend.position="none",
          axis.title = element_text(size =12),
          legend.text = element_text(size =12),
          legend.title = element_text(size =12),
          strip.text.x = element_text(size = 12),
          strip.background = element_rect(fill="white"))
 box}

cpg1<-plt_of_CpG_val_less("cg07839457", "TI")
cpg1
ggsave(file=here("DNAm/figs","NLRC5_rescopes_boxonly.pdf"),cpg1, width=5, height=3.5)
ggsave(file=here("DNAm/figs/jpeg","NLRC5_rescopes_boxonly.jpeg"),cpg1, width=5, height=3.5)


plt_of_CpG_val_less_cor<-function(CpG, tissue){
  smp_info_plt<-sampleinfo_DNAm[which(sampleinfo_DNAm$sample.site==tissue),]
  smp_info_plt$betas<-combat_ibd_Beta[CpG,which(sampleinfo_DNAm$sample.site==tissue)]
  smp_info_plt$diagnosis_grouped<-factor(smp_info_plt$diagnosis_grouped, levels=c("Control","UC","CD"))
  sampleinfo_re_ID<-sampleinfo_DNAm[which(!(is.na(sampleinfo_DNAm$passage.or.rescope.no))),"sample_ID"]
  sampleinfo_DNAm_reog<-sampleinfo_DNAm[which(sampleinfo_DNAm$sample_ID%in%sampleinfo_re_ID),]
  ibd_combo_reog<-combat_ibd_Beta[,which(colnames(combat_ibd_Beta)%in%sampleinfo_DNAm_reog$array.id)]
  smp_info_plt<-sampleinfo_DNAm_reog[which(sampleinfo_DNAm_reog$sample.site==tissue),]
  smp_info_plt$betas<-ibd_combo_reog[CpG,which(sampleinfo_DNAm_reog$sample.site==tissue)]
  og<-smp_info_plt[which(is.na(smp_info_plt$passage.or.rescope.no)),]
  re<-smp_info_plt[which(!(is.na(smp_info_plt$passage.or.rescope.no))),]
  re<-re[which(re$case.no%in%og$case.no),]
  re<-re[order(re$case.no),]
  og<-og[order(og$case.no),]
  plt<-data.frame(case.no=re$case.no, re_beta=re$betas, og_beta=og$betas)
  plt$mean<-rowMeans(plt[,2:3])
  plt$diff<-plt$og_beta-plt$re_beta
  
  corcpg=signif(cor.test(re$betas,og$betas, method="spearman")$estimate, 2)
  pvalcpg=signif(cor.test(re$betas,og$betas, method="spearman")$p.value, 2)
  
  plt<-merge(plt, smp_info_plt[,c(2,21)], by="case.no")
  
  cor<-ggplot(plt, aes(og_beta, re_beta))+geom_point(aes(fill=diagnosis_grouped),shape=21, color="black", size=2)+fillscale_diagnosis+
    stat_smooth(method="lm", se=F, color="black")+ylim(0,1)+xlim(0,1)+theme_bw()+annotate("text", x=0.5, y=0.8, label= paste("Rs = ",corcpg,";","p = ",pvalcpg))+
    theme(legend.position = "none")+xlab("Original Sample DNAm")+ylab("Rescope Sample DNAm")
 }

cpg1<-plt_of_CpG_val_less_cor("cg07839457", "TI")
cpg1
ggsave(file=here("DNAm/figs","NLRC5_rescopes_cornly.pdf"),cpg1, width=3, height=2.9)
ggsave(file=here("DNAm/figs/jpeg","NLRC5_rescopes_coronly.jpeg"),cpg1, width=3, height=2.9)








#' ## stats only function

# CpGs<-TI_sig_UCCTRL$CpG
# tissue<-"TI"
# 
# CpGs<-SC_CD$CpG
# tissue<-"SC"

stats_function<-function(CpGs, tissue){
  do.call(rbind, lapply(1:length(CpGs), function(x){
    CpG<-CpGs[x]
    sampleinfo_re_ID<-sampleinfo_DNAm[which(!(is.na(sampleinfo_DNAm$passage.or.rescope.no))),"sample_ID"]
    sampleinfo_DNAm_reog<-sampleinfo_DNAm[which(sampleinfo_DNAm$sample_ID%in%sampleinfo_re_ID),]
    ibd_combo_reog<-combat_ibd_Beta[,which(colnames(combat_ibd_Beta)%in%sampleinfo_DNAm_reog$array.id)]
    identical(colnames(ibd_combo_reog),sampleinfo_DNAm_reog$array.id)
    smp_info_plt<-sampleinfo_DNAm_reog[which(sampleinfo_DNAm_reog$sample.site==tissue),]
    smp_info_plt$betas<-ibd_combo_reog[CpG,which(sampleinfo_DNAm_reog$sample.site==tissue)]
    og<-smp_info_plt[which(is.na(smp_info_plt$passage.or.rescope.no)),]
    re<-smp_info_plt[which(!(is.na(smp_info_plt$passage.or.rescope.no))),]
    re<-re[which(re$case.no%in%og$case.no),]
    re<-re[order(re$case.no),]
    og<-og[order(og$case.no),]
    identical(re$case.no, og$case.no)
    
    data.frame(CpG=CpG,cor=cor(re$betas,og$betas, method="spearman"), rmse=rmse(re$betas,og$betas))
  }))
  }

SC_CD_rescope<-stats_function(SC_CD$CpG, "SC")
mean(SC_CD_rescope$cor)
mean(SC_CD_rescope$rmse)

TI_CD_rescope<-stats_function(TI_sig_UCCTRL$CpG, "TI")
mean(TI_CD_rescope$cor)
mean(TI_CD_rescope$rmse)

SC_UC_rescope<-stats_function(SC_UC$CpG, "SC")
mean(SC_UC_rescope$cor)
mean(SC_UC_rescope$rmse)



##################
#'# Differential Methylation Analysis
##################
## rescopes and controls
sampleinfo_re<-sampleinfo_DNAm[which(!(is.na(sampleinfo_DNAm$passage.or.rescope.no))),]
sampleinfo_ctrl<-sampleinfo_DNAm[which(sampleinfo_DNAm$diagnosis_grouped=="Control"),]

sampleinfo_rectrl<-rbind(sampleinfo_re,sampleinfo_ctrl)


#subset tissue
sampleinfo_rectrl_TI<-sampleinfo_rectrl[which(sampleinfo_rectrl$sample.site=="TI"),]
ibd_Mval_rectrl_TI<-ibd_Mval[,which(colnames(ibd_Mval)%in%sampleinfo_rectrl_TI$array.id)]
sampleinfo_rectrl_TI<-sampleinfo_rectrl_TI[match(colnames(ibd_Mval_rectrl_TI),sampleinfo_rectrl_TI$array.id),]
ibd_combo_rectrl_TI<-combat_ibd_Beta[,which(colnames(combat_ibd_Beta)%in%sampleinfo_rectrl_TI$array.id)]
identical(colnames(ibd_Mval_rectrl_TI),sampleinfo_rectrl_TI$array.id)
identical(colnames(ibd_combo_rectrl_TI),sampleinfo_rectrl_TI$array.id)

sampleinfo_rectrl_SC<-sampleinfo_rectrl[which(sampleinfo_rectrl$sample.site=="SC"),]
ibd_Mval_rectrl_SC<-ibd_Mval[,which(colnames(ibd_Mval)%in%sampleinfo_rectrl_SC$array.id)]
sampleinfo_rectrl_SC<-sampleinfo_rectrl_SC[match(colnames(ibd_Mval_rectrl_SC),sampleinfo_rectrl_SC$array.id),]
ibd_combo_rectrl_SC<-combat_ibd_Beta[,which(colnames(combat_ibd_Beta)%in%sampleinfo_rectrl_SC$array.id)]
identical(colnames(ibd_Mval_rectrl_SC),sampleinfo_rectrl_SC$array.id)
identical(colnames(ibd_combo_rectrl_SC),sampleinfo_rectrl_SC$array.id)







#' ### SC CD
mval<-ibd_Mval_rectrl_SC[which(rownames(ibd_Mval_rectrl_SC)%in%SC_CD$CpG),]
beta<-ibd_combo_rectrl_SC[which(rownames(ibd_combo_rectrl_SC)%in%SC_CD$CpG),]
meta<-sampleinfo_rectrl_SC

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
SC_CD_re<-limmares

SC_CD_compare<-merge(SC_CD_re[,c(1,2,7)],SC_CD[,c(1,2,7)], by="CpG")
colnames(SC_CD_compare)<-c("CpG","db_rescope","pval_rescope","db_original","pval_original")

ggplot(SC_CD_compare, aes(db_original, db_rescope))+geom_hline(yintercept=c(-0.05,0.05), color="grey")+geom_vline(xintercept=c(-0.05,0.05), color="grey")+
  geom_point(shape=21,size=2, color="black", fill="#a6d96a")+th_present+theme_bw()+xlab("Original Delta Beta")+ylab("Rescope Delta Beta")+ylim(-0.2,0.2)+xlim(-0.2,0.2)

ggsave(file=here("DNAm/figs","SC_CD_rescope_DB.pdf"), width=4,height=3.75)
ggsave(file=here("DNAm/figs/jpeg","SC_CD_rescope_DB.jpeg"), width=4,height=3.75)

cor(SC_CD_compare$db_rescope, SC_CD_compare$db_original, method="spearman")


#' How many significant? CD
## sig num strict
stat_hits<-as.data.frame(SC_CD_re)[which(SC_CD_re$adjusted_p_CD<=0.05),]
SC_CD_rescope<-stat_hits[which(abs(stat_hits$db_ctrl_CD)>=0.05),]  #14/52
nrow(SC_CD_rescope)

#' How many significant? UC
## sig num strict
SC_UC_compare<-merge(SC_CD_re[,c(1,3,6)],SC_CD[,c(1,3,6)], by="CpG")
colnames(SC_UC_compare)<-c("CpG","db_rescope","pval_rescope","db_original","pval_original")

stat_hits<-as.data.frame(SC_CD_re)[which(SC_CD_re$adjusted_p_CD<=0.05),]
SC_CD_rescope<-stat_hits[which(abs(stat_hits$db_ctrl_CD)>=0.05),]  #14/52
nrow(SC_CD_rescope)





#'### TI CD
mval<-ibd_Mval_rectrl_TI[which(rownames(ibd_Mval_rectrl_TI)%in%TI_sig_UCCTRL$CpG),]
beta<-ibd_combo_rectrl_TI[which(rownames(ibd_combo_rectrl_TI)%in%TI_sig_UCCTRL$CpG),]
meta<-sampleinfo_rectrl_TI

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

diff_meth_TI_UCCTRL_re <- data.frame(p.value=ebfit$p.value[,"diagnosis_TI_splitCD"], delta_beta=delta_beta)
diff_meth_TI_UCCTRL_re$CpG<-rownames(diff_meth_TI_UCCTRL_re)
diff_meth_TI_UCCTRL$CpG<-rownames(diff_meth_TI_UCCTRL)


TI_CD_compare<-merge(diff_meth_TI_UCCTRL_re,diff_meth_TI_UCCTRL, by="CpG")
colnames(TI_CD_compare)<-c("CpG","pval_rescope","db_rescope","pval_original","db_original")
ggplot(TI_CD_compare, aes(db_original, db_rescope))+geom_hline(yintercept=c(-0.05,0.05), color="grey")+geom_vline(xintercept=c(-0.05,0.05), color="grey")+
  geom_point(shape=21,size=2, color="black", fill="cornflowerblue")+th_present+theme_bw()+xlab("Original Delta Beta")+ylab("Rescope Delta Beta")+ylim(-0.2,0.2)+xlim(-0.2,0.2)

ggsave(file=here("DNAm/figs","TI_CD_rescope_DB.pdf"), width=4,height=3.75)
ggsave(file=here("DNAm/figs/jpeg","TI_CD_rescope_DB.jpeg"), width=4,height=3.75)

cor(TI_CD_compare$db_rescope, TI_CD_compare$db_original, method="spearman")

diff_meth_TI_UCCTRL_re$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL_re$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL_re))

### How many significant
# sig num strict
stat_hits<-as.data.frame(diff_meth_TI_UCCTRL_re)[which(diff_meth_TI_UCCTRL_re$adjusted_p<=0.05),]
TI_sig_UCCTRL_re<-stat_hits[which(abs(stat_hits$delta_beta)>=0.05),] # 61/110
nrow(TI_sig_UCCTRL_re)

# CpG hisglighted in paper and significant in primary (cg11706729 was on epic but only sig in organoid)
#btoh sig in rescopes
TI_sig_UCCTRL_re[grep("cg02756056|cg07839457",TI_sig_UCCTRL_re$CpG),]

## How many MHC are maintained longtidunally
MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")
EPIC_genes<-read.csv(here("DNAm/data","EPIC_ensembl_gene_annotation.csv")) # 1137194
#EPIC_genes<-read.csv(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","EPIC_ensembl_gene_annotation.csv")) # 1137194
MHC_CpG<-EPIC_genes[which(EPIC_genes$Gene.name%in%MHCI),]

TI_sig_UCCTRL_re_MHC<-TI_sig_UCCTRL_re[which(TI_sig_UCCTRL_re$CpG%in%MHC_CpG$IlmnID),]
nrow(TI_sig_UCCTRL_re_MHC)
diff_meth_TI_UCCTRL_re_MHC<-diff_meth_TI_UCCTRL_re[which(diff_meth_TI_UCCTRL_re$CpG%in%MHC_CpG$IlmnID),]
nrow(diff_meth_TI_UCCTRL_re_MHC)

              
              ### plot db highlight MHC I
              TI_CD_compare$MHCI<-"Not MHC I"
              TI_CD_compare[which(TI_CD_compare$CpG%in%MHC_CpG$IlmnID),]$MHCI<-"MHC I"
              ggplot(TI_CD_compare, aes(db_original, db_rescope, color=MHCI))+geom_hline(yintercept=c(-0.05,0.05), color="grey")+geom_vline(xintercept=c(-0.05,0.05), color="grey")+
                geom_point(shape=21,size=2,  fill="cornflowerblue")+th_present+theme_bw()+xlab("Original Delta Beta")+ylab("Rescope Delta Beta")+ylim(-0.2,0.2)+xlim(-0.2,0.2)+
                scale_color_manual(values=c("black","white"))
              ggsave(file=here("DNAm/figs","TI_CD_rescope_DB_MHCI_highlight.pdf"), width=4.5,height=3.25)
              ggsave(file=here("DNAm/figs/jpeg","TI_CD_rescope_DB_MHCI_highlight.jpeg"), width=4.5,height=3.25)
              
              
              SC_CD_compare$MHCI<-"Not MHC I"
              #SC_CD_compare[which(SC_CD_compare$CpG%in%MHC_CpG$IlmnID),]$MHCI<-"MHC I" # no MHC I in SC CD IEC
              ggplot(SC_CD_compare, aes(db_original, db_rescope, color=MHCI))+geom_hline(yintercept=c(-0.05,0.05), color="grey")+geom_vline(xintercept=c(-0.05,0.05), color="grey")+
                geom_point(shape=21,size=2, fill="#a6d96a")+th_present+theme_bw()+xlab("Original Delta Beta")+ylab("Rescope Delta Beta")+ylim(-0.2,0.2)+xlim(-0.2,0.2)+
                scale_color_manual(values=c("white"))
              ggsave(file=here("DNAm/figs","SC_CD_rescope_DB_MHCI_highlight.pdf"), width=4.5,height=3.25)
              ggsave(file=here("DNAm/figs/jpeg","SC_CD_rescope_DB_MHCI_highlight.jpeg"), width=4.5,height=3.25)
              
              cor(SC_CD_compare$db_rescope, SC_CD_compare$db_original, method="spearman")
              

#' ### SC UC
mval<-ibd_Mval_rectrl_SC[which(rownames(ibd_Mval_rectrl_SC)%in%SC_UC$CpG),]
beta<-ibd_combo_rectrl_SC[which(rownames(ibd_combo_rectrl_SC)%in%SC_UC$CpG),]
meta<-sampleinfo_rectrl_SC

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
avebeta.lm<-lm(beta~as.factor(inflammation)+as.factor(sex)+age, data=meta)
residuals<-residuals(summary(avebeta.lm))
colnames(residuals)<-colnames(beta)
adj.residuals<-residuals+matrix(mean(beta), nrow=1, ncol=length(residuals))

delta_beta<-lapply(1:nrow(adj.residuals), function(x) {
  group_means<-tapply(adj.residuals[x,], meta$diagnosis_grouped, mean)
  c(group_means[1]-group_means[2],group_means[1]-group_means[3], group_means[3]-group_means[2])
})

limmares<-as.data.frame(do.call(rbind, delta_beta))
colnames(limmares)<-c("db_ctrl_CD","db_ctrl_UC","db_UC_CD")
limmares$CpG<-SC_UC$CpG

p.value_UC<-data.frame(p.value_UC=contrastCpG_UC$P.Value, CpG=SC_UC$CpG)
p.value_CD<-data.frame(p.value_CD=contrastCpG_CD$P.Value, CpG=SC_UC$CpG)
p.value_CDUC<-data.frame(p.value_CDUC=contrastCpG_CDUC$P.Value, CpG=SC_UC$CpG)
pval<-merge(p.value_CDUC,p.value_UC, by="CpG")
pval<-merge(pval,p.value_CD, by="CpG")
limmares<-merge(limmares,pval, by="CpG")

limmares$adjusted_p_CDUC<-p.adjust(limmares$p.value_CDUC, method = "fdr", n = nrow(limmares))
limmares$adjusted_p_UC<-p.adjust(limmares$p.value_UC, method = "fdr", n = nrow(limmares))
limmares$adjusted_p_CD<-p.adjust(limmares$p.value_CD, method = "fdr", n = nrow(limmares))
SC_UC_re<-limmares

SC_UC_compare<-merge(SC_UC_re[,c(1,3,6)],SC_UC[,c(1,3,6)], by="CpG")
colnames(SC_UC_compare)<-c("CpG","db_rescope","pval_rescope","db_original","pval_original")

#' Not going to plot one point


#' How many significant?
## sig num strict
stat_hits<-as.data.frame(SC_UC_re)[which(SC_UC_re$adjusted_p_UC<=0.05),]
SC_UC_rescope<-stat_hits[which(abs(stat_hits$db_ctrl_UC)>=0.05),]  #14/52
nrow(SC_UC_rescope)




#'## R Session Info
sessionInfo()




# 
# ################ ################ ################ ################ ################ ################ ################ TI ALL CpG
# ##### TI
# mval<-ibd_Mval_rectrl_TI
# beta<-ibd_combo_rectrl_TI
# meta<-sampleinfo_rectrl_TI
# 
# meta$diagnosis_grouped<-as.factor(meta$diagnosis_grouped)
# meta$diagnosis_TI_split<-meta$diagnosis_grouped
# levels(meta$diagnosis_TI_split)<-c("Control_UC","CD","Control_UC")
# 
# meta$inflammation<-as.factor(meta$inflammation)
# meta$sex<-as.factor(meta$sex)
# 
# mod<-model.matrix(~diagnosis_TI_split+inflammation+sex+age, data=meta)
# fit <- lmFit(mval, mod)
# ebfit <- eBayes(fit)
# 
# avebeta.lm<-apply(beta, 1, function(x){
#   lm(x~as.factor(inflammation)+as.factor(sex)+age, data=meta)})
# residuals<-t(sapply(avebeta.lm, function(x)residuals(summary(x))))
# colnames(residuals)<-colnames(beta)
# adj.residuals<-residuals+matrix(apply(beta, 1, mean), nrow=nrow(residuals), ncol=ncol(residuals))
# 
# delta_beta<-sapply(1:nrow(adj.residuals), function(x) {
#   group_means<-tapply(adj.residuals[x,], meta$diagnosis_TI_split, mean)
#   group_means[1]-group_means[2] # Control_UC - CD
# })
# 
# diff_meth_TI_UCCTRL_re <- data.frame(p.value=ebfit$p.value[,"diagnosis_TI_splitCD"], delta_beta=delta_beta)
# diff_meth_TI_UCCTRL_re$CpG<-rownames(diff_meth_TI_UCCTRL_re)
# diff_meth_TI_UCCTRL$CpG<-rownames(diff_meth_TI_UCCTRL)
# 
# 
# TI_CD_compare<-merge(diff_meth_TI_UCCTRL_re,diff_meth_TI_UCCTRL, by="CpG")
# colnames(TI_CD_compare)<-c("CpG","pval_rescope","db_rescope","pval_original","db_original")
# save(TI_CD_compare, file="../../output/TI_rescope_all_stats.RData")
# 
# load( file="../../output/TI_rescope_all_stats.RData")
# ggplot(TI_CD_compare, aes(db_original, db_rescope))+geom_hline(yintercept=c(-0.05,0.05), color="grey")+geom_vline(xintercept=c(-0.05,0.05), color="grey")+
#   geom_point(shape=21,size=2, color="black", fill="cornflowerblue")+th+theme_bw()+xlab("Delta Beta Original")+ylab("Delta Beta Rescope")+ylim(-0.2,0.2)+xlim(-0.2,0.2)
# 
# 
# TI_CD_compare$adjusted_p_rescope<-p.adjust(TI_CD_compare$pval_rescope, method = "fdr", n = nrow(TI_CD_compare))
# TI_CD_compare$adjusted_p_original<-p.adjust(TI_CD_compare$pval_original, method = "fdr", n = nrow(TI_CD_compare))
# 
# 
# print(cor(TI_CD_compare$db_rescope, TI_CD_compare$db_original, method="spearman"))
# 
# TI_CD_compare$color<-sapply(1:nrow(TI_CD_compare), function(x){
#   if(TI_CD_compare$adjusted_p_original[x] < 0.05 & abs(TI_CD_compare$db_original[x])>0.05 & TI_CD_compare$adjusted_p_rescope[x] < 0.05 & abs(TI_CD_compare$db_rescope[x])>0.05){"Original and Rescope"}else{
#     if(TI_CD_compare$adjusted_p_original[x] < 0.05 & abs(TI_CD_compare$db_original[x])>0.05){"Original"}else{
#       if(TI_CD_compare$adjusted_p_rescope[x] < 0.05 & abs(TI_CD_compare$db_rescope[x])>0.05){"Rescope"}else{"Not\nSignificant"}
#     }
#   }
# })
# 
# TI_CD_compare<-TI_CD_compare[(order(TI_CD_compare$color)),]
# 
# CpG_OI<-c("cg24512973","cg18804777")
# 
# TI_CD_compare$hightlight<-" "
# TI_CD_compare$hightlight[which(TI_CD_compare$CpG%in%CpG_OI)]<-"MUC1 CpG"
# 
# ggplot(TI_CD_compare, aes(db_original,db_rescope))+
#   geom_point(aes(fill=color, alpha=color, color=hightlight), shape=21,size=2)+theme_bw()+th+
#   geom_hline(yintercept=0, color="grey30")+
#   geom_vline(xintercept=0, color="grey30")+
#   geom_vline(xintercept=c(0.05, -0.05), color="grey70")+
#   geom_hline(yintercept=c(0.05, -0.05), color="grey70")+
#   xlim(-0.20, 0.2)+ylim(-0.2, 0.2)+
#   scale_fill_manual(values=c("lightgrey","#99ccff","#1874CD","#65b3ff"), name="Differential\nDNAm")+
#   scale_alpha_manual(values=c(0.25,1,1,1), name="Differential\nDNAm")+
#   scale_color_manual(values=c("white","black"), name="")+
#   stat_smooth(method="lm", color="black")+
#   xlab("Original Delta Beta")+ ylab("Rescope Delta Beta")+th
# 
# ggsave("../../figs/jpeg/delta_beta_cor_TICD_rescope_original_MUC1.jpeg", width = 7, height = 4.5)
