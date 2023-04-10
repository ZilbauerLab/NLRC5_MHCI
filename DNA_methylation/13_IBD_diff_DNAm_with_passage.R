#'---
#'title: Organoid and Primary differenital DNAm comparison
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
#library(limma)
library(gridExtra)
library(here)

options(stringsAsFactors = FALSE)


#' functions
source(here("general_functions","00_pretty_plots.R"))


manhattan_data<-function(CpG, pvalue, fdr, db, diagnosis, tissue, sample.type){
  data.frame(CpG=CpG, db=db, p.value=pvalue, fdr=fdr, diagnosis=diagnosis, tissue=tissue, sample.type=sample.type)
}
sig_hits<-function(manhattan_data_df){
  manhattan_data_df[which(manhattan_data_df$fdr<0.05 & abs(manhattan_data_df$db)>0.05),]
}



#' ## Load DNAm data
load(file=paste(here("DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))
# load(file=paste(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))

## fix 202 and 203
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605870061_R03C01")]<-"hold"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605880096_R04C01")]<-"205605870061_R03C01"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="hold")]<-"205605880096_R04C01"




#' Plot heteroskedacsity over passage

plt_hetero<-function(CpGs, legend, axislab, title){
  betas<-melt(combat_organoid_Beta[CpGs,])
  epic.organoid_plt<-merge(epic.organoid_combined, betas, by.x="array.id",by.y="X2")
  
  p<-ggplot(epic.organoid_plt, aes(passage.or.rescope.no_numeric,value))+
    geom_line(aes(group=sample_ID),color="lightgrey")+
    stat_smooth(method="lm", color="grey30", size=0.7, se=F)+th+theme_bw()+
    geom_point(aes(fill=as.factor(passage.or.rescope.no_numeric)),shape=21, size=1.25)+
    scale_fill_brewer(palette = "Spectral",name="Passage\nNumber")+facet_wrap(~X1)+
    ylab("DNAm Beta")+xlab("Passage Number")+ylim(0,1)+
    theme(plot.margin = margin(0.5, 0.15, 0.5, 0.15, "cm"),plot.title = element_text(size=12))
  
  if(missing(legend) & missing(axislab) & missing(title)){p}else{
    if(legend=="N" & axislab=="N"){p + theme(legend.position = "none",axis.title.y=element_blank(),
                                             axis.text.y=element_blank(),
                                             axis.ticks.y=element_blank())+ ggtitle(title)}else{
                                               if(legend=="N" & axislab=="Y"){p + theme(legend.position = "none") + ggtitle(title)}}}
}

plt_hetero(c("cg07839457","cg07862320"))
ggsave(here("DNAm/figs/jpeg","NLRC5_Passage.jpeg"), width = 6.5, height = 3.5)
ggsave(here("DNAm/figs","NLRC5_Passage.pdf"), width = 6.5, height = 3.5)





#' ## Load differential CpG from Organoids
load(here("DNAm/data/organoid_diff_DNAm_grouped_casectrl_ThreeBatch.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/organoid_diff_DNAm_grouped_casectrl_ThreeBatch.RData"))
diff_meth_TI_organoid<-diff_meth_TI_grouped
diff_meth_SC_organoid<-diff_meth_SC_grouped

manhattan_TI_UC_organoid<-manhattan_data(diff_meth_TI_organoid$CpG, diff_meth_TI_organoid$p.value_UC, diff_meth_TI_organoid$adjusted_p_UC, diff_meth_TI_organoid$db_ctrl_UC, "UC", "TI", "organoid")
manhattan_SC_UC_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_UC, diff_meth_SC_organoid$adjusted_p_UC, diff_meth_SC_organoid$db_ctrl_UC, "UC", "SC", "organoid")
manhattan_SC_CD_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_CD, diff_meth_SC_organoid$adjusted_p_CD, diff_meth_SC_organoid$db_ctrl_CD, "CD", "SC", "organoid")

load(here("DNAm/data/diff_DNAm_TI_UCCTRL_organoid_ThreeBatch.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/diff_DNAm_TI_UCCTRL_organoid_ThreeBatch.RData"))

diff_meth_TI_UCCTRL_organoid$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL_organoid$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL_organoid))

manhattan_TI_CD_organoid<-manhattan_data(rownames(diff_meth_TI_UCCTRL_organoid), diff_meth_TI_UCCTRL_organoid$p.value, diff_meth_TI_UCCTRL_organoid$adjusted_p, diff_meth_TI_UCCTRL_organoid$delta_beta, "CD", "TI", "organoid")


#In SC compare IBD vs ctrl
load(here("DNAm/data","organoid_diff_DNAm_simple_casectrl_ThreeBatch.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","organoid_diff_DNAm_simple_casectrl_ThreeBatch.RData"))

manhattan_SC_UCCD_organoid<-manhattan_data(rownames(diff_meth_SC_organoid), diff_meth_SC_organoid$p.value, diff_meth_SC_organoid$adjusted_p, diff_meth_SC_organoid$delta_beta, "IBD", "SC", "organoid")



###############
### load update diagnosis
###############
epic.samples<-read.table(here("DNAm/data","AllEPICOrganoid_MHCIannotated_UpdatedSampleInfo_31Oct22.txt"), header=T, sep="\t")
#epic.samples<-read.table(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","AllEPICOrganoid_MHCIannotated_UpdatedSampleInfo_31Oct22.txt"), header=T, sep="\t")

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



#'## plot effect of passage
SC_UC<-sig_hits(manhattan_SC_UC_organoid)
SC_CD<-sig_hits(manhattan_SC_CD_organoid)
TI_CD<-sig_hits(manhattan_TI_CD_organoid)
SC_IBD<-sig_hits(manhattan_SC_UCCD_organoid)


plt_passage_diff<-function(CpG, tissue){
  meta<-epic.organoid_exclude_other
  meta_tissue<-meta[which(meta$sample.site==tissue),]
  beta_tissue<-combat_organoid_Beta[,which(meta$sample.site==tissue)]
  
  meta_tissue$beta<-beta_tissue[CpG,]
  
  meta_tissue$diagnosis<-factor(meta_tissue$diagnosis, levels=c("Control", "UC","CD"))
  
  
  myColors_diagnosis_pass <- c("lightgrey","#FFB90F","dodgerblue3",
                               "lightgrey","#FFB90F","dodgerblue3",
                               "lightgrey","#FFB90F","dodgerblue3",
                               "lightgrey","#FFB90F","dodgerblue3",
                               "lightgrey","#FFB90F","dodgerblue3",
                               "lightgrey","#FFB90F","dodgerblue3",
                               "lightgrey","#FFB90F","dodgerblue3",
                               "lightgrey","#FFB90F","dodgerblue3",
                               "lightgrey","#FFB90F","dodgerblue3")
                               
  color_possibilities_diagnosis_pass<- levels(interaction(meta_tissue$diagnosis, meta_tissue$passage.or.rescope.no_numeric))
  names(myColors_diagnosis_pass) <- color_possibilities_diagnosis_pass
  fillscale_diagnosis_pass <- scale_fill_manual(name="Diagnosis and Passage",
                                           values = myColors_diagnosis_pass, drop = F, guide=F)
  txt<-meta_tissue[,c(8, 15)]
  txt<-txt[!duplicated(txt),]
  
  meta_tissue<-meta_tissue[which(meta_tissue$diagnosis%in%c("Control", "UC","CD")),]
  
  ggplot(meta_tissue, aes(diagnosis, beta, fill=interaction(diagnosis, passage.or.rescope.no_numeric)))+
    geom_boxplot(outlier.shape=NA)+geom_point(position = position_dodge(width=0.75), shape=21, size=1,color="black")+
    th+theme_bw()+ylim(0,1)+ylab("DNA Methylation")+xlab("Diagnosis")+fillscale_diagnosis_pass+
    geom_text(aes(label=passage.or.rescope.no_numeric, y=1),data=txt,position = position_dodge(width=0.75), color="grey60")
  }

#' TI_CD hit (gene:HS1BP3)
plt_passage_diff("cg17605007","TI")
ggsave(file=here("DNAm/figs/jpeg","CpG_passage_effect_cg17605007.jpeg"), plt_passage_diff("cg17605007","TI"))
ggsave(file=here("DNAm/figs","CpG_passage_effect_cg17605007.pdf"), plt_passage_diff("cg17605007","TI"))

#' TI_CD hit (gene:RAD51B)
plt_passage_diff("cg11469373","TI")
ggsave(file=here("DNAm/figs/jpeg","CpG_passage_effect_cg11469373.jpeg"), plt_passage_diff("cg11469373","TI"))
ggsave(file=here("DNAm/figs","CpG_passage_effect_cg11469373.pdf"), plt_passage_diff("cg11469373","TI"))


#' NLRC5
grid.arrange(plt_passage_diff("cg07839457","TI"),plt_passage_diff("cg07862320","TI"))
ggsave(file=here("DNAm/figs/jpeg","CpG_passage_effect_NLRC5.jpeg"), grid.arrange(plt_passage_diff("cg07839457","TI"),plt_passage_diff("cg07862320","TI")))
ggsave(file=here("DNAm/figs","CpG_passage_effect_NLRC5.pdf"), grid.arrange(plt_passage_diff("cg07839457","TI"),plt_passage_diff("cg07862320","TI")))




#######################
#'# passage and organoid differential CpG overlap
#######################

#' from previous analysis (just copied over)
load(file=here("DNAm/data","Heteroskedactiy_pvalues_FDR_1000iter_w_CpG.Rdata"))
load(file=here("DNAm/data","Heteroskedactiy_pvalues_FDR_1000iter_w_CpG.Rdata"))
#load(file="/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/Heteroskedactiy_pvalues_FDR_1000iter_w_CpG.Rdata")

pvals_long$BP_pval<-((1000-pvals_long$BP_count)+1)/1001
pvals_long$diff_pval<-((1000-pvals_long$diff_count)+1)/1001
pvals_long$BP_fdr<-((1000-pvals_long$fdr_BP)+1)/1001
pvals_long$diff_fdr<-((1000-pvals_long$fdr_diff)+1)/1001

pvals_long<-pvals_long[which(pvals_long$CpG%in%rownames(combat_organoid_Beta)),]
ibd_EPIC_organoid_beta_passage<-combat_organoid_Beta[which(rownames(combat_organoid_Beta)%in%pvals_long$CpG),]

pvals_long<-pvals_long[match(pvals_long$CpG,rownames(ibd_EPIC_organoid_beta_passage)),]
identical(rownames(ibd_EPIC_organoid_beta_passage), pvals_long$CpG)

length(diff_CpG<-rownames(ibd_EPIC_organoid_beta_passage)[which(pvals_long$diff_fdr<0.05)])#27341
diff_CpG_db<-pvals_long[which(pvals_long$diff_fdr<0.05 & abs(pvals_long$mean_db)>0.15),] #23761

length(intersect(SC_CD$CpG, diff_CpG_db$CpG)) #524/7410
length(intersect(SC_UC$CpG, diff_CpG_db$CpG)) #7/711
length(intersect(TI_CD$CpG, diff_CpG_db$CpG)) #26/522
length(intersect(SC_IBD$CpG, diff_CpG_db$CpG)) #506/7883


#' ### Genes effected by passage
EPIC_genes<-read.csv(here("DNAm/data","EPIC_ensembl_gene_annotation.csv")) # 1137194
#EPIC_genes<-read.csv(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","EPIC_ensembl_gene_annotation.csv")) # 1137194


unique(EPIC_genes[which(EPIC_genes$IlmnID%in%intersect(SC_CD$CpG, diff_CpG_db$CpG)),"Gene.name"])
unique(EPIC_genes[which(EPIC_genes$IlmnID%in%intersect(SC_UC$CpG, diff_CpG_db$CpG)),"Gene.name"])
unique(EPIC_genes[which(EPIC_genes$IlmnID%in%intersect(TI_CD$CpG, diff_CpG_db$CpG)),"Gene.name"])
unique(EPIC_genes[which(EPIC_genes$IlmnID%in%intersect(SC_IBD$CpG, diff_CpG_db$CpG)),"Gene.name"])

#' ### How many MHCI
MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")


sig_CpG<-rbind(SC_CD, SC_UC, TI_CD, SC_IBD)
IBD_genes<-merge(sig_CpG, EPIC_genes[,c(3,6,7,9)], by.x="CpG", by.y="IlmnID")

IBD_genes[which(IBD_genes$Gene.name%in%MHCI & IBD_genes$CpG%in%diff_CpG_db$CpG),]

MHCI[which(MHCI%in%unique(EPIC_genes[which(EPIC_genes$IlmnID%in%intersect(TI_CD$CpG, diff_CpG_db$CpG)),"Gene.name"]))]
MHCI[which(MHCI%in%unique(EPIC_genes[which(EPIC_genes$IlmnID%in%intersect(SC_IBD$CpG, diff_CpG_db$CpG)),"Gene.name"]))]

#'## R Session Info
sessionInfo()


