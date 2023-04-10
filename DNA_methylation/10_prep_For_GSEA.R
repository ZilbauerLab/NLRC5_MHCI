#'---
#'title: Prep lists fro Go ORA
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---
#'

#' #### Load Libraries
library(here)
library(dplyr)
suppressMessages(library(reshape))
library(ggplot2)
library(tidytext)

options(stringsAsFactors = FALSE)

source(here("general_functions","00_pretty_plots.R"))


#' ## functions
manhattan_data<-function(CpG, pvalue, fdr, db, diagnosis, tissue, sample.type){
data.frame(CpG=CpG, db=db, p.value=pvalue, fdr=fdr, diagnosis=diagnosis, tissue=tissue, sample.type=sample.type)
}

sig_hits<-function(manhattan_data_df){
  manhattan_data_df[which(manhattan_data_df$fdr<0.05 & abs(manhattan_data_df$db)>0.05),]
}

#'## Load CpG stats

load(file=here("DNAm/data/diff_DNAm_grouped_casectrl.RData"))
#load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/diff_DNAm_grouped_casectrl.RData"))

diff_meth_TI_primary<-diff_meth_TI_grouped
diff_meth_SC_primary<-diff_meth_SC_grouped

manhattan_TI_UC_primary<-manhattan_data(diff_meth_TI_primary$CpG, diff_meth_TI_primary$p.value_UC,diff_meth_TI_primary$adjusted_p_UC, diff_meth_TI_primary$db_ctrl_UC, "UC", "TI", "primary")
manhattan_SC_UC_primary<-manhattan_data(diff_meth_SC_primary$CpG, diff_meth_SC_primary$p.value_UC,diff_meth_SC_primary$adjusted_p_UC, diff_meth_SC_primary$db_ctrl_UC, "UC", "SC", "primary")
manhattan_SC_CD_primary<-manhattan_data(diff_meth_SC_primary$CpG, diff_meth_SC_primary$p.value_CD,diff_meth_SC_primary$adjusted_p_CD, diff_meth_SC_primary$db_ctrl_CD, "CD", "SC", "primary")

load(here("DNAm/data/organoid_diff_DNAm_grouped_casectrl_ThreeBatch.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/organoid_diff_DNAm_grouped_casectrl_ThreeBatch.RData"))

diff_meth_TI_organoid<-diff_meth_TI_grouped
diff_meth_SC_organoid<-diff_meth_SC_grouped

manhattan_TI_UC_organoid<-manhattan_data(diff_meth_TI_organoid$CpG, diff_meth_TI_organoid$p.value_UC, diff_meth_TI_organoid$adjusted_p_UC, diff_meth_TI_organoid$db_ctrl_UC, "UC", "TI", "organoid")
manhattan_SC_UC_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_UC, diff_meth_SC_organoid$adjusted_p_UC, diff_meth_SC_organoid$db_ctrl_UC, "UC", "SC", "organoid")
manhattan_SC_CD_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_CD, diff_meth_SC_organoid$adjusted_p_CD, diff_meth_SC_organoid$db_ctrl_CD, "CD", "SC", "organoid")

load(file=here("DNAm/data/diff_DNAm_TI_UCCTRL.RData"))
#load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/diff_DNAm_TI_UCCTRL.RData"))

load(here("DNAm/data/diff_DNAm_TI_UCCTRL_organoid_ThreeBatch.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/diff_DNAm_TI_UCCTRL_organoid_ThreeBatch.RData"))

diff_meth_TI_UCCTRL$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL))
diff_meth_TI_UCCTRL_organoid$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL_organoid$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL_organoid))

manhattan_TI_CD_primary<-manhattan_data(rownames(diff_meth_TI_UCCTRL), diff_meth_TI_UCCTRL$p.value, diff_meth_TI_UCCTRL$adjusted_p, diff_meth_TI_UCCTRL$delta_beta, "CD", "TI", "primary")
manhattan_TI_CD_organoid<-manhattan_data(rownames(diff_meth_TI_UCCTRL_organoid), diff_meth_TI_UCCTRL_organoid$p.value, diff_meth_TI_UCCTRL_organoid$adjusted_p, diff_meth_TI_UCCTRL_organoid$delta_beta, "CD", "TI", "organoid")


#In SC compare IBD vs ctrl
load(here("DNAm/data","organoid_diff_DNAm_simple_casectrl_ThreeBatch.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","organoid_diff_DNAm_simple_casectrl_ThreeBatch.RData"))
manhattan_SC_UCCD_organoid<-manhattan_data(rownames(diff_meth_SC_organoid), diff_meth_SC_organoid$p.value, diff_meth_SC_organoid$adjusted_p, diff_meth_SC_organoid$delta_beta, "IBD", "SC", "organoid")

load(here("DNAm/data","diff_DNAm_simple_casectrl.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","diff_DNAm_simple_casectrl.RData"))
manhattan_SC_UCCD_primary<-manhattan_data(rownames(diff_meth_SC), diff_meth_SC$p.value, diff_meth_SC$adjusted_p, diff_meth_SC$delta_beta, "IBD", "SC", "primary")


## EPIC only CpGs
load(file=here("DNAm/data","diff_DNAm_EPIC_only.RData"))
#load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","diff_DNAm_EPIC_only.RData"))

diff_meth_TI_UCCTRL$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL))

manhattan_TI_CD_primary_EPIC<-manhattan_data(rownames(diff_meth_TI_UCCTRL), diff_meth_TI_UCCTRL$p.value, diff_meth_TI_UCCTRL$adjusted_p, diff_meth_TI_UCCTRL$delta_beta, "CD", "TI", "primary")
manhattan_SC_UC_primary_EPIC<-manhattan_data(diff_meth_SC_grouped$CpG, diff_meth_SC_grouped$p.value_UC,diff_meth_SC_grouped$adjusted_p_UC, diff_meth_SC_grouped$db_ctrl_UC, "UC", "SC", "primary")
manhattan_SC_CD_primary_EPIC<-manhattan_data(diff_meth_SC_grouped$CpG, diff_meth_SC_grouped$p.value_CD,diff_meth_SC_grouped$adjusted_p_CD, diff_meth_SC_grouped$db_ctrl_CD, "CD", "SC", "primary")
manhattan_SC_UCCD_primary_EPIC<-manhattan_data(rownames(diff_meth_SC), diff_meth_SC$p.value,diff_meth_SC$adjusted_p, diff_meth_SC$delta_beta, "IBD", "SC", "primary")


rm(diff_meth_AC_grouped, diff_meth_SC_grouped, diff_meth_TI_grouped, diff_meth_TI_UCCTRL)




#' ## Gene CpG annotation
# Made in previous scripts
EPIC_genes<-read.csv(here("DNAm/data","EPIC_ensembl_gene_annotation.csv")) # 1137194
#EPIC_genes<-read.csv(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","EPIC_ensembl_gene_annotation.csv")) # 1137194



#########
## pathway enrichment
#########
source("DNAm/scripts/GSEA_function_tmp.R")
GO_file = here("DNAm/data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

plot_GSEA<-function(CpG_diff_DNAm){
  CpG_diff_DNAm_genes<-merge(CpG_diff_DNAm, EPIC_genes, by.x="CpG", by.y="IlmnID")
  CpG_diff_DNAm_genes<-CpG_diff_DNAm_genes[,c("Gene.name","db")]
  CpG_diff_DNAm_genes_max<-as.data.frame(CpG_diff_DNAm_genes %>% group_by(Gene.name) %>% summarise_each(funs(.[which.max(abs(.))])))
  
  gene_list = CpG_diff_DNAm_genes_max$db
  names(gene_list) = CpG_diff_DNAm_genes_max$Gene.name
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  
  res = GSEA(gene_list, GO_file, pval = 0.05)
  
  plt_path<-res$Results
  plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
  plt_path$Enrichment_Cell<-"Up-regulated in \nrecently recruited myeloid"
  plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nKupffer-like cells"
  
  plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
  
  plt_path$direction_label<-as.factor(plt_path$Enrichment)
  levels(plt_path$direction_label)<-c(0.1,-0.1)
  plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
  
  # top 10
  plt_path[1:10,]}



TI_CD_GSEA<-plot_GSEA(manhattan_TI_CD_organoid)
TI_CD_GSEA$comparison<-"Hypomethylated in CD compared to Controls + UC in TI"
SC_CD_GSEA<-plot_GSEA(manhattan_SC_CD_organoid)
SC_CD_GSEA$comparison<-"Hypomethylated in CD compared to Controls in SC"
SC_UC_GSEA<-plot_GSEA(manhattan_SC_UC_organoid)
SC_UC_GSEA$comparison<-"Hypomethylated in UC compared to Controls in SC"


plt_path<-rbind(TI_CD_GSEA, SC_UC_GSEA)




GSEA_CD_UC<-plt_path %>%
  group_by(comparison) %>%
  mutate(comparison = as.factor(comparison),
         pathway = reorder_within(pathway, NES, comparison)) %>%
  ggplot(aes(NES, pathway))+geom_point(aes(size=size, fill=comparison), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+scale_y_reordered()+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_vline(xintercept = 0, color="white")+scale_fill_manual(values=c("dodgerblue3","darkgoldenrod1"))+
  guides(fill = F,size = guide_legend(title = "Geneset\nSize"))+
  facet_wrap(~comparison, ncol=1, scales="free_y")+
  theme(legend.position = c(0.12, 0.84),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="grey40"),
        panel.grid = element_blank(),strip.background = element_rect(colour=NA, fill=NA))

GSEA_CD_UC
ggsave(file="DNAm/figs/GSEA_TICD_UCSC.pdf",GSEA_CD_UC,w=12,h=6)




plt_path<-rbind(TI_CD_GSEA, SC_CD_GSEA ,SC_UC_GSEA)

GSEA_CD_UC_SCTI<-plt_path %>%
  group_by(comparison) %>%
  mutate(comparison = as.factor(comparison),
         pathway = reorder_within(pathway, NES, comparison)) %>%
  ggplot(aes(NES, pathway))+geom_point(aes(size=size, fill=comparison), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+scale_y_reordered()+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_vline(xintercept = 0, color="white")+scale_fill_manual(values=c("dodgerblue3","dodgerblue3","darkgoldenrod1"))+
  guides(fill = F,size = guide_legend(title = "Geneset\nSize"))+
  facet_wrap(~comparison, ncol=1, scales="free_y")+
  theme(legend.position = c(0.12, 0.84),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="grey40"),
        panel.grid = element_blank(),strip.background = element_rect(colour=NA, fill=NA))

GSEA_CD_UC_SCTI
ggsave(file="DNAm/figs/GSEA_TICD_SCCD_UCSC.pdf",GSEA_CD_UC_SCTI,w=12,h=9)








MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","MR1","CD1D","IRF1","NLRC5")

#########
## MHC I focused pathway enrichment
#########
source("DNAm/scripts/GSEA_function_tmp.R")
GO_file = here("DNAm/data/Human_GOBP_AllPathways_with_GO_iea_October_26_2022_symbol.gmt")

plot_GSEA_MHCI<-function(CpG_diff_DNAm){
  sig_hits_list<-sig_hits(CpG_diff_DNAm)
  
  CpG_diff_DNAm_genes<-merge(CpG_diff_DNAm, EPIC_genes, by.x="CpG", by.y="IlmnID")
  sig_genes<-CpG_diff_DNAm_genes[which(CpG_diff_DNAm_genes$CpG%in%sig_hits_list$CpG),]
  sig_genes<-sig_genes[which(sig_genes$db>0),]
  sig_genes<-unique(sig_genes$Gene.name)
  
  CpG_diff_DNAm_genes<-CpG_diff_DNAm_genes[,c("Gene.name","db")]
  CpG_diff_DNAm_genes_max<-as.data.frame(CpG_diff_DNAm_genes %>% group_by(Gene.name) %>% summarise_each(funs(.[which.max(abs(.))])))
  
  gene_list = CpG_diff_DNAm_genes_max$db
  names(gene_list) = CpG_diff_DNAm_genes_max$Gene.name
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  
  res = GSEA(gene_list, GO_file, pval = 0.05)
  
  plt_path<-res$Results
  plt_path$pathway<-sapply(1:nrow(plt_path), function(x) strsplit(plt_path$pathway[x], "%")[[1]][1])
  plt_path$Enrichment_Cell<-"Up-regulated in \nrecently recruited myeloid"
  plt_path$Enrichment_Cell[which(plt_path$Enrichment=="Down-regulated")]<-"Up-regulated in \nKupffer-like cells"
  
  #plt_path$label<-lapply(1:nrow(plt_path), function(x) paste0(plt_path$leadingEdge[x][[1]][1:4], collapse = ", "))
  plt_path$label<-lapply(1:nrow(plt_path), function(x) {
    gene_label<-plt_path$leadingEdge[x][[1]]
    sig_gene_label<-gene_label[which(gene_label%in%sig_genes)]
    not_sig_gene_label<-gene_label[which(!(gene_label%in%sig_genes))]
    MHCI_gene_label<-sig_gene_label[which(sig_gene_label%in%MHCI)]
    not_MHCI_gene_label<-sig_gene_label[which(!(sig_gene_label%in%MHCI))]
    if(length(MHCI_gene_label)>=4){paste0(MHCI_gene_label[1:4], collapse = ", ")}else{
      if(length(MHCI_gene_label)==0){
        if(length(sig_gene_label)<4){paste0(c(sig_gene_label, not_sig_gene_label[1:(4-length(sig_gene_label))]), collapse = ", ")}else{
          paste0(sig_gene_label[1:4], collapse = ", ")}}else{
            if(length(MHCI_gene_label)<4){paste0(c(MHCI_gene_label, not_MHCI_gene_label[1:(4-length(MHCI_gene_label))]), collapse = ", ")}
          }
    }
  })
  
  plt_path$direction_label<-as.factor(plt_path$Enrichment)
  levels(plt_path$direction_label)<-c(0.1,-0.1)
  plt_path$direction_label<-as.numeric(as.character(plt_path$direction_label))
  
  # top 10
  plt_path[1:10,]}




TI_CD_GSEA<-plot_GSEA_MHCI(manhattan_TI_CD_organoid)
TI_CD_GSEA$comparison<-"Hypomethylated in CD compared to Controls + UC in TI"
SC_CD_GSEA<-plot_GSEA_MHCI(manhattan_SC_CD_organoid)
SC_CD_GSEA$comparison<-"Hypomethylated in CD compared to Controls in SC"
SC_UC_GSEA<-plot_GSEA_MHCI(manhattan_SC_UC_organoid)
SC_UC_GSEA$comparison<-"Hypomethylated in UC compared to Controls in SC"


plt_path<-rbind(TI_CD_GSEA, SC_CD_GSEA ,SC_UC_GSEA)

GSEA_CD_UC_SCTI<-plt_path %>%
  group_by(comparison) %>%
  mutate(comparison = as.factor(comparison),
         pathway = reorder_within(pathway, NES, comparison)) %>%
  ggplot(aes(NES, pathway))+geom_point(aes(size=size, fill=comparison), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+scale_y_reordered()+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_vline(xintercept = 0, color="white")+scale_fill_manual(values=c("dodgerblue3","dodgerblue3","darkgoldenrod1"))+
  guides(fill = F,size = guide_legend(title = "Geneset\nSize"))+
  facet_wrap(~comparison, ncol=1, scales="free_y")+
  theme(legend.position = c(0.12, 0.84),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="grey40"),
        panel.grid = element_blank(),strip.background = element_rect(colour=NA, fill=NA))

GSEA_CD_UC_SCTI
ggsave(file="DNAm/figs/GSEA_TICD_SCCD_UCSC_MHCI_highlight.pdf",GSEA_CD_UC_SCTI,w=12,h=9)




########### overlap GSEA
overlap_manhattan<-merge(manhattan_TI_CD_organoid, manhattan_SC_CD_organoid, by="CpG")
overlap_manhattan$db<-sapply(1:nrow(overlap_manhattan), function(x) mean(overlap_manhattan$db.x[x],overlap_manhattan$db.y[x]))
overlap_manhattan$fdr<-sapply(1:nrow(overlap_manhattan), function(x) max(overlap_manhattan$fdr.x[x],overlap_manhattan$fdr.y[x]))
overlap_manhattan<-overlap_manhattan[,c("CpG","db","fdr")]

overlap_CD_GSEA<-plot_GSEA_MHCI(overlap_manhattan)
overlap_CD_GSEA$comparison<-"Hypomethylated in CD compared to Controls in both TI and SC"

plt_path<-rbind(overlap_CD_GSEA ,SC_UC_GSEA)

GSEA_CD_UC_SCTI<-plt_path %>%
  group_by(comparison) %>%
  mutate(comparison = as.factor(comparison),
         pathway = reorder_within(pathway, NES, comparison)) %>%
  ggplot(aes(NES, pathway))+geom_point(aes(size=size, fill=comparison), shape=21)+
  theme_bw()+th_present+ylab("")+xlab("Normalized Enrichment Score")+scale_y_reordered()+
  geom_text(aes(label=label),hjust="inward",  nudge_x = plt_path$direction_label, color="grey50", size=3)+
  geom_vline(xintercept = 0, color="white")+scale_fill_manual(values=c("dodgerblue3","darkgoldenrod1"))+
  guides(fill = F,size = guide_legend(title = "Geneset\nSize"))+
  facet_wrap(~comparison, ncol=1, scales="free_y")+
  theme(legend.position = c(0.12, 0.84),
        legend.background = element_rect(size=0.5, linetype="solid", 
                                         colour ="grey40"),
        panel.grid = element_blank(),strip.background = element_rect(colour=NA, fill=NA))

GSEA_CD_UC_SCTI
ggsave(file="DNAm/figs/GSEA_TISCCD_UCSC_MHCI_highlight.pdf",GSEA_CD_UC_SCTI,w=12,h=6)
