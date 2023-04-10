
#'### Load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(here)
library(ggplot2)
library(reshape2)
library(gridExtra)
#library(limma)
library(cowplot)
library(RColorBrewer)
library(ggsignif)


load(here("../EBI/scRNAseq_codon/data/kong_immunity","Kong_Immunity_epithelial_MHCI_raw.RData"))

TI_UMAP<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/cluster/TI_EPI.scp.X_umap.coords.txt", header=T)
TI_UMAP<-TI_UMAP[-1,]
CO_UMAP<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/cluster/CO_EPI.scp.X_umap.coords.txt", header=T)
CO_UMAP<-CO_UMAP[-1,]

meta<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/metadata/scp_metadata_combined.v2.txt",header=T)
meta<-meta[-1,]



##############
### TI
##############
plt<-merge(meta, TI_UMAP, by="NAME")
plt<-merge(plt, score_data_allcompartment, by="NAME")

plt<-plt[which(plt$Site == "TI"),]
mn_sample_MHCI<-as.data.frame(tapply(plt$MHCI_score1_allcompartment, plt$donor_id, mean))
colnames(mn_sample_MHCI)<-"Mean_MHCI_score"
mn_sample_MHCI$individual<-rownames(mn_sample_MHCI)
mn_sample_MHCI$individual<-as.factor(mn_sample_MHCI$individual)

meta_min<-meta[,c("donor_id","disease__ontology_label")]
meta_min<-meta_min[!duplicated(meta_min),]

meta_TI<-meta[which(meta$Site=="TI"),]
cell_type_count<-table(meta_TI$donor_id, meta_TI$Celltype)
cell_type_percent<-(cell_type_count/rowSums(cell_type_count))*100

cell_type_percent<-as.data.frame(cell_type_percent)


cell_count_score<-merge(cell_type_percent, mn_sample_MHCI, by.x="Var1",by.y="individual")
cell_count_score<-merge(cell_count_score, meta_min, by.x="Var1",by.y="donor_id")


summary_stat<-as.data.frame(cell_count_score[grep("T cells",cell_count_score$Var2),] %>%  
                              group_by(Var2, disease__ontology_label) %>% summarise(cor_coef = cor.test(Mean_MHCI_score, Freq, method="spearman")$estimate,
                                                                       p_val = cor.test(Mean_MHCI_score, Freq, method="spearman")$p.value))

ggplot(cell_count_score, aes(Mean_MHCI_score, Freq))+geom_point()+facet_wrap(~Var2, scales="free_y")

ggplot(cell_count_score[grep("T cells",cell_count_score$Var2),], 
       aes(Mean_MHCI_score, Freq, fill=disease__ontology_label, color=disease__ontology_label))+
  geom_point(color="black",shape=21)+facet_grid(disease__ontology_label~Var2, scales="free_y")+
  stat_smooth(method="lm",se=F)+theme_bw()+theme(legend.position = "none")+
  geom_text(aes(x=0.4, y=15, label=paste("Rs = ",signif(cor_coef,2),";","p = ",signif(p_val,1))),data=summary_stat, color="black")
ggsave(file="../EBI/MHCI/MHCI_paper_figs/Tcell_proportion_MHCIscore_KongTI.pdf", w=20,h=10)


##############
### SC
##############
plt<-merge(meta, CO_UMAP, by="NAME")
plt<-merge(plt, score_data_allcompartment, by="NAME")

mn_sample_MHCI<-as.data.frame(tapply(plt$MHCI_score1_allcompartment, plt$donor_id, mean))
colnames(mn_sample_MHCI)<-"Mean_MHCI_score"
mn_sample_MHCI$individual<-rownames(mn_sample_MHCI)
mn_sample_MHCI$individual<-as.factor(mn_sample_MHCI$individual)

meta_min<-meta[,c("donor_id","disease__ontology_label")]
meta_min<-meta_min[!duplicated(meta_min),]

meta_CO<-meta[which(meta$Site=="CO"),]
cell_type_count<-table(meta_CO$donor_id, meta_CO$Celltype)
cell_type_percent<-(cell_type_count/rowSums(cell_type_count))*100

cell_type_percent<-as.data.frame(cell_type_percent)


cell_count_score<-merge(cell_type_percent, mn_sample_MHCI, by.x="Var1",by.y="individual")
cell_count_score<-merge(cell_count_score, meta_min, by.x="Var1",by.y="donor_id")


summary_stat<-as.data.frame(cell_count_score[grep("T cells",cell_count_score$Var2),] %>%  
                              group_by(Var2, disease__ontology_label) %>% summarise(cor_coef = cor.test(Mean_MHCI_score, Freq, method="spearman")$estimate,
                                                                                    p_val = cor.test(Mean_MHCI_score, Freq, method="spearman")$p.value))

ggplot(cell_count_score, aes(Mean_MHCI_score, Freq))+geom_point()+facet_wrap(~Var2, scales="free_y")

ggplot(cell_count_score[grep("T cells",cell_count_score$Var2),], 
       aes(Mean_MHCI_score, Freq, fill=disease__ontology_label, color=disease__ontology_label))+
  geom_point(color="black",shape=21)+facet_grid(disease__ontology_label~Var2, scales="free_y")+
  stat_smooth(method="lm",se=F)+theme_bw()+theme(legend.position = "none")+
  geom_text(aes(x=0.4, y=15, label=paste("Rs = ",signif(cor_coef,2),";","p = ",signif(p_val,1))),data=summary_stat, color="black")
ggsave(file="../EBI/MHCI/MHCI_paper_figs/Tcell_proportion_MHCIscore_KongSC.pdf", w=20,h=10)

