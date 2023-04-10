#'---
#'title: Smilie sc data analysis
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

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

options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))


#' #### color pallet
#' ### cell labels
cols_manual<-c( "#00441b","#6521C1",
                "#238b45", "#7fbc41","#41ab5d","#74c476",
                "#B921C1","#4eb3d3","#308AC8", "#C86E30","#C12134",
                "#fec44f","#fe9929","#ec7014","#993404",
                "#053061")

names(cols_manual) <- c("Best4+ Enterocytes","Stem",
                        "Enterocytes","Enterocyte Progenitors","Immature Enterocytes 2", "Immature Enterocytes 1",
                        "Enteroendocrine","Immature Goblet", "Goblet" , "Paneth cell","Tuft",
                        "TA 1","TA 2","Cycling TA","Secretory TA",
                        "M cells")

fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)



#' ## processed data 
#' Downloaded here:https://singlecell.broadinstitute.org/single_cell/study/SCP259/intra-and-inter-cellular-rewiring-of-the-human-colon-during-ulcerative-colitis#/
#' Raw is available with request

dataset_loc <- here("data/smilie/SCP259/expression")

d10x.epi <- Read10X(file.path("data/smilie/SCP259/expression/epi"), gene.column=1)
#' Initialize the Seurat object with the raw (non-normalized data).
d10x.epi <- CreateSeuratObject(counts = d10x.epi, min.cells = 3, min.features = 0)


meta<-read.csv(here("data/smilie/SCP259/metadata/all.meta2.txt"), sep="\t", header=T)
meta<-meta[-1,]

#meta_cell_add<-meta_cell_add[match(meta_cell_add$cell, colnames(d10x)),]
meta_cell_add<-meta[match(colnames(d10x.epi), meta$NAME),]
rm(meta)
identical(meta_cell_add$NAME, rownames(d10x.epi@meta.data))
rownames(meta_cell_add)<-meta_cell_add$NAME

d10x.epi<- AddMetaData(d10x.epi, meta_cell_add)


# meta from supplment
meta_diagnosis<-read.csv(here("data/smilie/supp_S1.csv"), header=T)
#Don't have scRNA seq data for N897 (so only 17 UC when 18 reported), but have data for N110 and no information on diagnosis from paper supplement
# do have all 12 healthy controls
meta_diagnosis<-rbind(meta_diagnosis,data.frame(Subject.ID="N110",Disease="Unknown",Location="Unknown",Gender="Unknown",Smoking="Unknown"))
meta_diagnosis_add<-meta_diagnosis[match(d10x.epi@meta.data$Subject, meta_diagnosis$Subject.ID),]

identical(meta_diagnosis_add$Subject.ID, d10x.epi@meta.data$Subject)

rownames(meta_diagnosis_add)<-d10x.epi@meta.data$NAME

d10x.epi<- AddMetaData(d10x.epi, meta_diagnosis_add)
rm(meta_diagnosis_add)


########################
#' Gene scores
########################

MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")
crypt_villis = c("SEPP1", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8", "APOBEC1", "APOB", "APOA4", "APOA1", "NPC1L1", "EGFR", "KLF4", "ENPP3", "NT5E", "SLC28A2", "ADA")


##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.epi <- NormalizeData(d10x.epi,scale.factor = 10000, normalization.method = "LogNormalize")

d10x.epi <- AddModuleScore(
  object = d10x.epi,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score'
)

d10x.epi <- AddModuleScore(
  object = d10x.epi,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)

plt<-d10x.epi@meta.data[,c("MHCI_score1","crypt_villis_score1","Cluster","Disease","Health")]

plt<-plt[which(plt$Disease!="Unknown"),]






# Group the data by cell type do another anova
stats_plt_stim<-lapply(unique(plt$Cluster), function(x){
  print("##################")
  cell_meta<-plt[which(plt$Cluster==x),]
  print(unique(as.character(cell_meta$Cluster)))
  pairwise<-pairwise.t.test(cell_meta$MHCI_score1, cell_meta$Disease, pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$orig.ident<-x
  
  mns<-tapply(cell_meta$MHCI_score1, cell_meta$Disease, mean)
  
  pairwise$mn_diff<-mns[1]-mns[2]#UC-healthy
  pairwise
})

stats_plt_stim<-do.call(rbind, stats_plt_stim)
stats_plt_stim$Var1<-as.character(stats_plt_stim$Var1)
stats_plt_stim$Var2<-as.character(stats_plt_stim$Var2)

stats_plt_stim$p_adjusted<-p.adjust(stats_plt_stim$value, method="fdr", n=nrow(stats_plt_stim))
sig_UC_ctrl<-stats_plt_stim[which(stats_plt_stim$p_adjusted<0.05),]
sig_UC_ctrl


mn_cryptvillus<-tapply(plt$crypt_villis_score1, plt$Cluster, mean)
cryptvillus_order<-names(mn_cryptvillus)[(order(mn_cryptvillus))]
plt$Cluster<-factor(plt$Cluster, (cryptvillus_order))

sig_df<-data.frame(Cluster=cryptvillus_order, sig="")
sig_df$sig[which(sig_df$Cluster%in%sig_UC_ctrl$orig.ident)]<-"*"
sig_df$direction<-"Lower UC"
sig_df$direction[which(sig_df$Cluster%in%sig_UC_ctrl[which(sig_UC_ctrl$mn_diff>0),"orig.ident"])]<-"Higher UC"

sig_df$Cluster<-factor(sig_df$Cluster, (cryptvillus_order))

ggplot(plt, aes(Disease,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.2,aes(fill=Disease))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~Cluster, nrow=2)+scale_fill_manual(values=c("darkgoldenrod1","lightgrey"))+scale_color_manual(values=c("darkgoldenrod1","grey40"), guide=F)+
  geom_text(data=sig_df, aes(x=1.5, y=1.75, label=sig, color=direction), size=8)

ggsave(file=here("figs", "smilie_MHCI_sortedby_crypt_villus.pdf"), w=15,h=6)
ggsave(file=here("figs/jpeg", "smilie_MHCI_sortedby_crypt_villus.jpeg"), w=15,h=6)




ggplot(plt, aes(Disease, MHCI_score1, fill=Disease))+geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.25)+
  scale_fill_manual(values=c("darkgoldenrod1","lightgrey"))+theme_bw()+th

t.test(plt$MHCI_score1~plt$Disease)
mns<-tapply(plt$MHCI_score1, plt$Disease, mean)
mns[1]-mns[2]#UC-healthy

ggsave(file=here("figs", "smilie_MHCI_all_epi.pdf"), w=6,h=6)
ggsave(file=here("figs/jpeg", "smilie_MHCI_all_epi.jpeg"), w=6,h=6)


########
#' # Crypt villus
########
plt$Cluster<-factor(plt$Cluster, levels = c("Stem", "Enterocytes","Enterocyte Progenitors","Immature Enterocytes 2", "Immature Enterocytes 1","Best4+ Enterocytes",
                                            "Enteroendocrine","Immature Goblet", "Goblet" , "Paneth cell","Tuft",
                                            "TA 1","TA 2","Cycling TA","Secretory TA","M cells"))
ggplot(plt, aes(reorder(Cluster, crypt_villis_score1),crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=Cluster))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"),
                      axis.text.x = element_text(angle = 45, hjust = 1))+
  fillscale_cols_manual+xlab("")

ggsave(file=here("figs", "smilie_crypt_villus_celltype.pdf"), w=15,h=6)
ggsave(file=here("figs/jpeg", "smilie_crypt_villus_celltype.jpeg"), w=15,h=6)



###### summary plot
scatter<-ggplot(plt, aes(crypt_villis_score1, MHCI_score1, fill=Disease, color=Disease))+geom_point(shape=21, color="black", size=1.5)+
  theme_bw()+th+theme(legend.position = "none",plot.margin = unit(c(0,0,1,4.5), "cm"))+
  scale_color_manual(values=c("darkgoldenrod1","lightgrey"))+scale_fill_manual(values=c("darkgoldenrod1","lightgrey"))+  stat_smooth(method="lm",se=F)


violin_celltype<-ggplot(plt, aes(reorder(Cluster, crypt_villis_score1), crypt_villis_score1, fill=Cluster))+
  geom_violin()+scale_x_discrete(position = "bottom") +coord_flip()+
  scale_fill_manual(values=cols_manual)+theme_bw()+th+ylab("")+xlab("")+theme(axis.title.x=element_blank(),
                                                                                   axis.text.x=element_blank(),
                                                                                   axis.ticks.x=element_blank(),
                                                                                   legend.position = "none",
                                                                                   plot.margin = unit(c(1,0,0,1), "cm"))

box_diagnosis<-ggplot(plt, aes(Disease, MHCI_score1, fill=Disease))+geom_boxplot()+
  scale_fill_manual(values=c("darkgoldenrod1","lightgrey"))+theme_bw()+th+theme(axis.title.y =element_blank(),
                                          axis.text.y=element_blank(),
                                          axis.ticks.y=element_blank(),
                                          plot.margin = unit(c(0,0.2,1,0), "cm"))

mn_cryptvillus<-tapply(plt$crypt_villis_score1, plt$Cluster, mean)
cryptvillus_order<-names(mn_cryptvillus)[(order(mn_cryptvillus))]

plt_count<-melt(table(plt$Disease, plt$Cluster))
plt_count<-plt_count[which(plt_count$Var2!="Paneth cell"),]
colnames(plt_count)<-c("Diagnosis","Cluster","count")
plt_count$cluster_ID<-factor(plt_count$Cluster, cryptvillus_order)
placeholder<-ggplot(plt_count, aes(Diagnosis, Cluster))+geom_tile(color="grey50", fill="white",size=0.4)+
  geom_text(aes(label=count),size=2.75,color="grey40")+ggtitle(" Cell Number")+
  theme_void()+theme(axis.title=element_blank(),
                     axis.text=element_blank(),
                     axis.ticks=element_blank(),
                     plot.margin = unit(c(0.4,3.1,0,0), "cm"),
                     plot.title = element_text(size = 10))
placeholder

grid_fig<-plot_grid(violin_celltype, placeholder, scatter, box_diagnosis, 
                    #align = 'hv',
                    ncol = 2, axis="lr", 
                    rel_heights = c(1,2),
                    rel_widths = c(4,1))

ggsave(file=here("figs", "smilie_MHCI_crpyt_summary.pdf"),grid_fig, w=10,h=10)
ggsave(file=here("figs/jpeg", "smilie_MHCI_crpyt_summary.jpeg"),grid_fig, w=10,h=10)

cor(plt$MHCI_score1, plt$crypt_villis_score1, method = "spearman")


#' 
#' ################
#' #' ## normalize and scale for umap
#' ################
#' d10x.epi <- NormalizeData(d10x.epi)
#' d10x.epi <- FindVariableFeatures(d10x.epi, selection.method = "vst", nfeatures = 2000)
#' d10x.epi <- ScaleData(d10x.epi) 
#' 
#' # dimension reduction
#' d10x.epi <- RunPCA(d10x.epi, ndims.print = 1:10, nfeatures.print = 10)
#' d10x.epi <- RunUMAP(d10x.epi, dims = 1:30)
#' #d10x.epi <- RunTSNE(d10x.epi, dims = 1:30)
#' 
#' 
#' DimPlot(d10x.epi, reduction = "umap", group.by = "Cluster", pt.size=0.5)
#' DimPlot(d10x.epi, reduction = "umap", group.by = "Health", pt.size=0.5)
#' DimPlot(d10x.epi, reduction = "umap", group.by = "Subject", pt.size=0.5)
#' DimPlot(d10x.epi, reduction = "umap", group.by = "Disease", pt.size=0.5)
#' 
#' 
#' 
#' d10x.epi@meta.data$Cluster<-factor(d10x.epi@meta.data$Cluster, levels = c("Stem",
#'                                                                            "Enterocytes","Enterocyte Progenitors","Immature Enterocytes 2", "Immature Enterocytes 1","Best4+ Enterocytes",
#'                                                                            "Enteroendocrine","Immature Goblet", "Goblet" , "Paneth cell","Tuft",
#'                                                                            "TA 1","TA 2","Cycling TA","Secretory TA",
#'                                                                            "M cells"))
#' 
#' DimPlot(d10x.epi, reduction = "umap", group.by = "Cluster", pt.size=0.5)+colcale_cols_manual
#' 



