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
library(here)
library(scales)
library(cowplot)
library(dplyr)


options(stringsAsFactors = FALSE)


#' functions
source(here("general_functions","00_pretty_plots.R"))
source(here("general_functions","00_DNAm_volcano.R"))





#' ## delta beta correlation plot
manhattan_data<-function(CpG, pvalue, fdr, db, diagnosis, tissue, sample.type){
  data.frame(CpG=CpG, db=db, p.value=pvalue, fdr=fdr, diagnosis=diagnosis, tissue=tissue, sample.type=sample.type)
}

sig_hits<-function(manhattan_data_df){
  manhattan_data_df[which(manhattan_data_df$fdr<0.05 & abs(manhattan_data_df$db)>0.05),]
}


#load(here("DNAm/data/organoid_diff_DNAm_grouped_casectrl_ThreeBatch.RData"))
load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/organoid_diff_DNAm_grouped_casectrl_ThreeBatch.RData"))

diff_meth_TI_organoid<-diff_meth_TI_grouped
diff_meth_SC_organoid<-diff_meth_SC_grouped

manhattan_TI_UC_organoid<-manhattan_data(diff_meth_TI_organoid$CpG, diff_meth_TI_organoid$p.value_UC, diff_meth_TI_organoid$adjusted_p_UC, diff_meth_TI_organoid$db_ctrl_UC, "UC", "TI", "organoid")
manhattan_TI_CD_organoid<-manhattan_data(diff_meth_TI_organoid$CpG, diff_meth_TI_organoid$p.value_CD, diff_meth_TI_organoid$adjusted_p_CD, diff_meth_TI_organoid$db_ctrl_CD, "CD", "TI", "organoid")
manhattan_TI_CDUC_organoid<-manhattan_data(diff_meth_TI_organoid$CpG, diff_meth_TI_organoid$p.value_CDUC, diff_meth_TI_organoid$adjusted_p_CDUC, diff_meth_TI_organoid$db_UC_CD, "CDUC", "TI", "organoid")

manhattan_SC_UC_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_UC, diff_meth_SC_organoid$adjusted_p_UC, diff_meth_SC_organoid$db_ctrl_UC, "UC", "SC", "organoid")
manhattan_SC_CD_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_CD, diff_meth_SC_organoid$adjusted_p_CD, diff_meth_SC_organoid$db_ctrl_CD, "CD", "SC", "organoid")
manhattan_SC_CDUC_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_CDUC, diff_meth_SC_organoid$adjusted_p_CDUC, diff_meth_SC_organoid$db_UC_CD, "CDUC", "SC", "organoid")


length(sig_hits(manhattan_TI_UC_organoid)$CpG)
length(sig_hits(manhattan_TI_CD_organoid)$CpG)
length(sig_hits(manhattan_TI_CDUC_organoid)$CpG)

intersect(sig_hits(manhattan_TI_UC_organoid)$CpG, sig_hits(manhattan_TI_CD_organoid)$CpG) #0
intersect(sig_hits(manhattan_TI_CDUC_organoid)$CpG, sig_hits(manhattan_TI_CD_organoid)$CpG) #0

length(intersect(sig_hits(manhattan_SC_UC_organoid)$CpG, sig_hits(manhattan_SC_CD_organoid)$CpG)) #577

length(sig_hits(manhattan_TI_CD_organoid)$CpG)#361
length(sig_hits(manhattan_SC_CD_organoid)$CpG)#7410
length(intersect(sig_hits(manhattan_TI_CD_organoid)$CpG, sig_hits(manhattan_SC_CD_organoid)$CpG)) #248


## how many hypomethylated
CD_overlap<-merge(sig_hits(manhattan_TI_CD_organoid), sig_hits(manhattan_SC_CD_organoid), by="CpG")
CD_overlap_hypo<-CD_overlap[which(CD_overlap$db.x>0 & CD_overlap$db.y>0),]

write.table(CD_overlap_hypo$CpG, file="DNAm/data/CD_hypomethylated_CpG_TI_SC.txt", row.names = F, quote = F, col.names = F)
