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



d10x.primary<-readRDS(here("data","d10x_primary_normalized.rds"))



## add cell type labels from split analysis
#immune
load(here("output","immune_iterative_label.Rdata"))
#stromal
load(here("output","strom_celltype_label.Rdata"))
#epithelial
load(here("output","epi_celltype_label.Rdata"))

cell_label<-rbind(epi_cell_labels, immune_cell_labels, stromal_cell_labels)

cell_label$cluster_ID<-as.character(cell_label$cluster_ID)

cell_label$cluster_ID[which(cell_label$cluster_ID=="Neonatal_cell")]<-"Neonatal Epithelial"
cell_label$cluster_ID[which(cell_label$cluster_ID=="crypt")]<-"Stem"


cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x.primary), cell_label$index),]
identical(colnames(d10x.primary), cell_label$index)


d10x.primary <- AddMetaData(d10x.primary, metadata = cell_label)

# ## rm unannotated
d10x.primary<-subset(d10x.primary, subset = cluster_ID != "unknown_filtered")

table(d10x.primary@meta.data$cluster_ID)


d10x.primary@meta.data$cell_type<-as.factor(d10x.primary@meta.data$cluster_ID)
d10x.primary@meta.data$cell_type<-factor(d10x.primary@meta.data$cell_type, levels = c(
                                                                                      "BEST4 enterocyte","Stem","early enterocyte","enterocyte","enteroendocrine",
                                                                                      "Goblet cell", "Paneth cell","Paneth (UC only)","Tuft",
  "Activated B cell","activated DC" ,"Activated T","Arterial endothelial cell" ,"B cell","CD4 T cell",
  "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
  "IgA plasma cell","IgG plasma cell","Lymphatic endothelial cell","Macrophage" ,"mast cells","Memory B cell" ,
  "Monocyte","pDC" ,"pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"TA","Tfh" ,
  "Treg", "Venous endothelial cell","Glial cell","myofibroblast","Mast cell"))

#' ### dev cell labels
cols_manual<-c(  
                 "#3D0AF2","#6521C1","#67D364","#367C34",
                 "#B921C1","#308AC8",
                 "#C86E30","#CFA218","#C12134","#238443","#810f7c" ,"#02818a","#8c510a" ,"#78c679","#67a9cf",
                 "#3690c0","#8073ac","#b2abd2","#238443","#dd3497","#1d91c0","#addd8e","#014636",
                 "#253494","#081d58","#bf812d","#7a0177" ,"#ce1256","#006d2c" ,
                 "#a50f15","#542788" ,"#e31a1c","#e08214","#ef6548","#fdd49e" ,"#f46d43","#4575b4" ,
                 "#016c59", "#bf812d","#c51b7d","#fb9a99","#fb6a4a")


names(cols_manual) <- levels(d10x.primary@meta.data$cell_type)
fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)



DimPlot(d10x.primary, reduction = "umap", group.by = "cell_type", pt.size=0.25)+colcale_cols_manual+ theme(legend.position="bottom")+ggtitle("")
  
ggsave(file=here("figs","primary_UMAP_celltype.pdf"), w=12,h=14)
ggsave(file=here("figs/jpeg","primary_UMAP_celltype.jpeg"), w=12,h=14)
