
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

source("/media/redgar/Seagate Portable Drive/EBI_backup/codon_july2022/redgar/scRNAseq_codon/scripts/00_pretty_plots.R")

MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")
crypt_villis = c("SEPP1", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8", "APOBEC1", "APOB", "APOA4", "APOA1", "NPC1L1", "EGFR", "KLF4", "ENPP3", "NT5E", "SLC28A2", "ADA")

## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
#d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))
d10x.primary<-readRDS("/media/redgar/Seagate Portable Drive/EBI_backup/codon_july2022/redgar/scRNAseq_codon/data/d10x_primary_raw_merged.rds")

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.primary <- NormalizeData(d10x.primary,scale.factor = 10000, normalization.method = "LogNormalize")

d10x.primary <- AddModuleScore(
  object = d10x.primary,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score'
)

d10x.primary <- AddModuleScore(
  object = d10x.primary,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)

score_data<-d10x.primary@meta.data[,c("MHCI_score1","crypt_villis_score1")]


# ## epithelial type labels
# load(here("data","primary.epi.cells.RData"))
load("/media/redgar/Seagate Portable Drive/EBI_backup/codon_july2022/redgar/scRNAseq_codon/data/primary.epi.cells.RData")

primary.epi.cells <- AddMetaData(primary.epi.cells, metadata = score_data)
#filter remaining immune (179 cells)
primary.epi.cells<-subset(primary.epi.cells, subset = cluster_ID != "B cell" & cluster_ID != "T cell"  & cluster_ID != "CD8 T cell")

d10x.primary_scores<-FetchData(object = d10x.primary, vars = c(MHCI))
d10x.primary_scores_epithelial<-d10x.primary_scores[which(rownames(d10x.primary_scores)%in%rownames(primary.epi.cells@meta.data)),]

d10x.primary_scores_epithelial$cell<-rownames(d10x.primary_scores_epithelial)
plt_counts<-melt(d10x.primary_scores_epithelial)

meta<-primary.epi.cells@meta.data
meta$cell<-rownames(meta)
plt_counts<-merge(plt_counts,meta[,c("cell","cluster_ID")] )

cols_manual_less<-c( "#87d435ff",
                     "cornflowerblue",
                     "#238b43ff","#5e1e6fff","#aa0044ff",
                     "#d9c026ff","#CFA218",  
                     "#f16916ff","#a13da1ff")

umap_mat<-as.data.frame(Embeddings(object = primary.epi.cells, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-primary.epi.cells@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")

plt_no_neo<-plt[which(plt$orig.ident!="Neonatal"),]
plt_no_neo<-plt_no_neo[which(plt_no_neo$cluster_ID!="Memory B cell"),]
plt_no_neo$cluster_ID<-as.character(plt_no_neo$cluster_ID)

plt_no_neo$cluster_ID<-as.factor(plt_no_neo$cluster_ID)
levels(plt_no_neo$cluster_ID)<-c("BEST4 Enterocyte","Stem","Enterocyte","Enteroendocrine","Goblet","Paneth","Paneth (UC only)","TA","Tuft")


mn_sample_MHCI<-as.data.frame(tapply(plt_no_neo$MHCI_score1, plt_no_neo$individual, mean))
colnames(mn_sample_MHCI)<-"Mean_MHCI_score"
mn_sample_MHCI$individual<-rownames(mn_sample_MHCI)
mn_sample_MHCI$individual<-as.factor(mn_sample_MHCI$individual)
levels(mn_sample_MHCI$individual)<-c("T017","T019","T024","T036","T036", "T036", "T057","T110", "T110", "T160",   
                          "T161","T175","T176","T180","T182","T184","T189","T197","T202","T203",   
                          "T44")

                    #load(here("/media/redgar/Seagate Portable Drive/EBI_backup/codon_july2022/redgar/scRNAseq_codon/output/cell_label_whole_clustering.RData"))
load("/home/redgar/Documents/EBI/scRNAseq_codon/data/Primary_UMAP_allsamples_plus_MHCI_raw.RData")
plt$cell_type<-as.character(plt$cell_type)
plt$cell_type[which(plt$cluster_ID=="Paneth")]<-"Paneth"
plt[which(is.na(plt$cell_type)),]

plt<-plt[which(plt$individual%in%c("T017","T019","T176","T189","T197","T202","T203","T024","T036","T44","T057","T160","T161","T175","T182","T184","T180")),]
plt$individual<-as.factor(plt$individual)


meta<-plt[,c("individual","orig.ident")]
meta<-meta[!duplicated(meta),]

# cell_type_count<-as.data.frame(tapply(plt$cell, list(plt$individual, plt$cell_type), length))
# cell_type_count[is.na(cell_type_count)]<-0
# cell_type_count$individual<-rownames(cell_type_count)
# cell_count_score<-merge(cell_type_count, mn_sample_MHCI, by="individual")


cell_type_count<-table(plt$individual, plt$cell_type)
cell_type_percent<-(cell_type_count/rowSums(cell_type_count))*100

cell_type_percent<-as.data.frame(cell_type_percent)


cell_count_score<-merge(cell_type_percent, mn_sample_MHCI, by.x="Var1",by.y="individual")
cell_count_score<-merge(cell_count_score, meta, by.x="Var1",by.y="individual")





summary_stat<-as.data.frame(cell_count_score[which(cell_count_score$Var2%in%c("Activated T","CD4 T cell","CD8 T cell") & cell_count_score$orig.ident!="UC"),] %>%  
  group_by(Var2, orig.ident) %>% summarise(cor_coef = cor.test(Mean_MHCI_score, Freq, method="spearman")$estimate,
                                                               p_val = cor.test(Mean_MHCI_score, Freq, method="spearman")$p.value))

ggplot(cell_count_score, aes(Mean_MHCI_score, Freq))+geom_point()+facet_wrap(~Var2, scales="free_y")

ggplot(cell_count_score[which(cell_count_score$Var2%in%c("Activated T","CD4 T cell","CD8 T cell")),], 
       aes(Mean_MHCI_score, Freq, fill=orig.ident, color=orig.ident))+
  geom_point(color="black",shape=21)+facet_grid(orig.ident~Var2, scales="free_y")+
  stat_smooth(method="lm",se=F)+
  fillscale_diagnosis+colscale_diagnosis+theme_bw()+th_present+theme(legend.position = "none")+
  geom_text(aes(x=0.4, y=25, label=paste("Rs = ",signif(cor_coef,2),";","p = ",signif(p_val,1))),data=summary_stat, color="black")
ggsave(file="../EBI/MHCI/MHCI_paper_figs/Tcell_proportion_MHCIscore.pdf", w=10,h=10)

summary_stat<-as.data.frame(cell_count_score[which(cell_count_score$Var2%in%c("Activated T","CD4 T cell","CD8 T cell")),] %>%  
                              group_by(Var2) %>% summarise(cor_coef = cor.test(Mean_MHCI_score, Freq, method="spearman")$estimate,
                                                                       p_val = cor.test(Mean_MHCI_score, Freq, method="spearman")$p.value))

ggplot(cell_count_score[which(cell_count_score$Var2%in%c("Activated T","CD4 T cell","CD8 T cell")),], 
       aes(Mean_MHCI_score, Freq))+
  geom_point(aes( fill=orig.ident),color="black",shape=21)+facet_wrap(~Var2, scales="free_y")+
  stat_smooth(method="lm",se=F, color="black")+
  fillscale_diagnosis+colscale_diagnosis+theme_bw()+th_present+theme(legend.position = "none")+
  geom_text(aes(x=0.3, y=10, label=paste("Rs = ",signif(cor_coef,2),";","p = ",signif(p_val,1))),data=summary_stat, color="black")




# 
# scatter<-ggplot(plt_no_neo, aes(crypt_villis_score1, MHCI_score1, fill=orig.ident, color=orig.ident))+geom_point(shape=21, color="black", size=1.5)+
#   fillscale_diagnosis+theme_bw()+th_present+theme(legend.position = "none",plot.margin = unit(c(0,0,1.95,3.65), "cm"))+scale_color_manual(values=c("dodgerblue3","grey50","darkgoldenrod3"))+
#   stat_smooth(method="lm",se=F)+xlab("Crypt-Villus Score")+ylab("MHC I Score")
# 
# 
# violin_celltype<-ggplot(plt_no_neo, aes(reorder(cluster_ID, crypt_villis_score1), crypt_villis_score1, fill=cluster_ID))+
#   geom_violin()+scale_x_discrete(position = "bottom") +coord_flip()+
#   scale_fill_manual(values=cols_manual_less)+
#   theme_bw()+th_present+ylab("")+xlab("")+theme(axis.title.x=element_blank(),
#                                                 axis.text.x=element_blank(),
#                                                 axis.ticks.x=element_blank(),
#                                                 legend.position = "none",
#                                                 plot.margin = unit(c(1,0,0,1), "cm"))
# 
# plt_no_neo$orig.ident<-as.character(plt_no_neo$orig.ident)
# table(plt_no_neo$orig.ident)
# box_diagnosis<-ggplot(plt_no_neo, aes(orig.ident, MHCI_score1, fill=orig.ident))+geom_boxplot()+xlab("Diagnosis")+
#   fillscale_diagnosis+theme_bw()+th_present+theme(axis.title.y =element_blank(),
#                                                   axis.text.y=element_blank(),
#                                                   axis.ticks.y=element_blank(),
#                                                   plot.margin = unit(c(0,0.2,1,0), "cm"),
#                                                   axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# mn_cryptvillus<-tapply(plt_no_neo$crypt_villis_score1, plt_no_neo$cluster_ID, mean)
# cryptvillus_order<-names(mn_cryptvillus)[(order(mn_cryptvillus))]
# 
# plt_count<-melt(table(plt_no_neo$orig.ident, plt_no_neo$cluster_ID))
# colnames(plt_count)<-c("orig.ident","cluster_ID","count")
# plt_count$cluster_ID<-factor(plt_count$cluster_ID, cryptvillus_order)
# placeholder<-ggplot(plt_count, aes(orig.ident, cluster_ID))+geom_tile(color="grey50", fill="white",size=0.4)+
#   geom_text(aes(label=count),size=3.75,color="grey40")+ggtitle(" Cell Number")+
#   theme_void()+theme(axis.title=element_blank(),
#                      axis.text=element_blank(),
#                      axis.ticks=element_blank(),
#                      plot.margin = unit(c(0.4,3.1,0,0), "cm"),
#                      plot.title = element_text(size = 12))
# placeholder
# 
# grid_fig<-plot_grid(violin_celltype, placeholder, scatter, box_diagnosis, 
#                     #align = 'hv',
#                     ncol = 2, axis="lr", 
#                     rel_heights = c(1,2),
#                     rel_widths = c(3.75,1.25))
# 
# ggsave(file=here("figs", "MHCI_crpyt_summary_present.pdf"),grid_fig, w=10,h=10)
# ggsave(file=here("figs/jpeg", "MHCI_crpyt_summary_present.jpeg"),grid_fig, w=10,h=10)
# 
