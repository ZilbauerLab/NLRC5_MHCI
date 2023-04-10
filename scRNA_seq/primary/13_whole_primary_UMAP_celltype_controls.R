#'---
#'title: scRNAseq seurat exploration
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
library(RColorBrewer)
library(ggsignif)


options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))


d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))

# relabel meta data
d10x.primary@meta.data$orig.ident[which(d10x.primary@meta.data$orig.ident=="Ctrl")]<-"Control"
d10x.primary@meta.data$orig.ident[which(d10x.primary@meta.data$orig.ident=="neo")]<-"Neonatal"

d10x.primary.control<-subset(d10x.primary, subset = orig.ident == "Control")


d10x.primary.control <- NormalizeData(d10x.primary.control)
d10x.primary.control <- FindVariableFeatures(d10x.primary.control, selection.method = "vst", nfeatures = 2000)
d10x.primary.control <- ScaleData(d10x.primary.control) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
d10x.primary.control <- RunPCA(d10x.primary.control, ndims.print = 1:10, nfeatures.print = 10)
d10x.primary.control <- RunUMAP(d10x.primary.control, dims = 1:30)
d10x.primary.control <- RunTSNE(d10x.primary.control, dims = 1:30)



######################
## cell cycle gene expression
######################

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x.primary.control <- CellCycleScoring(d10x.primary.control, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

## regress out cell cycle and other covariates
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#By default, sctransform accounts for cellular sequencing depth, or nUMIs.
d10x.primary.control <- SCTransform(d10x.primary.control, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE)


# dimension reduction
d10x.primary.control <- RunPCA(d10x.primary.control, verbose = FALSE)
d10x.primary.control <- RunUMAP(d10x.primary.control, dims = 1:30)
d10x.primary.control <- RunTSNE(d10x.primary.control, dims = 1:30)

# cluster
d10x.primary.control <- FindNeighbors(d10x.primary.control, reduction = "pca", dims = 1:20)
d10x.primary.control <- FindClusters(d10x.primary.control, resolution = 0.6)



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
cell_label<-cell_label[match(colnames(d10x.primary.control), cell_label$index),]
identical(colnames(d10x.primary.control), cell_label$index)


d10x.primary.control <- AddMetaData(d10x.primary.control, metadata = cell_label)

# ## rm unannotated
d10x.primary.control<-subset(d10x.primary.control, subset = cluster_ID != "unknown_filtered")

table(d10x.primary.control@meta.data$cluster_ID)


d10x.primary.control@meta.data$cell_type<-as.factor(d10x.primary.control@meta.data$cluster_ID)
d10x.primary.control@meta.data$cell_type<-factor(d10x.primary.control@meta.data$cell_type, levels = c(
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


names(cols_manual) <- levels(d10x.primary.control@meta.data$cell_type)
fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)

DimPlot(d10x.primary.control, reduction = "umap", group.by = "cell_type", pt.size=0.25)+colcale_cols_manual+ theme(legend.position="bottom")+ggtitle("")

ggsave(file=here("figs","primary_controls_UMAP_celltype.pdf"), w=4.5,h=6)
ggsave(file=here("figs/jpeg","primary_controls_UMAP_celltype.jpeg"), w=4.5,h=6)




umap_mat<-as.data.frame(Embeddings(object = d10x.primary.control, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.primary.control@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")

#########
## Controls only with compartment colors
#########
epithelial<-c("BEST4 enterocyte","Stem","early enterocyte","enterocyte","enteroendocrine",
              "Goblet cell", "Paneth cell","Paneth (UC only)","Tuft","TA","Paneth")
immune<-c("Activated B cell","activated DC" ,"Activated T","B cell","CD4 T cell",
          "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
          "IgA plasma cell","IgG plasma cell","Macrophage" ,"mast cells","Memory B cell" ,
          "Monocyte","pDC" ,"Tfh" ,
          "Treg","Mast cell")
stromal_and_glial<-c("S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"Lymphatic endothelial cell","Arterial endothelial cell","pericyte","Venous endothelial cell","myofibroblast", "Glial cell")

plt$general_type<-NA
plt$general_type[which(plt$cluster_ID%in%epithelial)]<-"Epithelial"
plt$general_type[which(plt$cluster_ID%in%immune)]<-"Immune"
plt$general_type[which(plt$cluster_ID%in%stromal_and_glial)]<-"Stromal and Glial"

getPalette = colorRampPalette(brewer.pal(9, "Reds"))
epi_col<-getPalette(length(unique(epithelial)))
getPalette = colorRampPalette(brewer.pal(9, "Blues"))
immune_col<-getPalette(length(unique(immune)))
getPalette = colorRampPalette(brewer.pal(9, "Greens"))
stromal_col<-getPalette(length(unique(stromal_and_glial)))

stromal_col[1]<-"#006D2C"
stromal_col[2]<-"#00441B"
stromal_col[8]<-"#41ae76"
stromal_col[9]<-"#66c2a4"
show_col(stromal_col)

epi_col[1]<-"#fc4e2a"
epi_col[2]<-"#a50f15"
show_col(epi_col)

immune_col[1]<-"#67a9cf"
immune_col[2]<-"#3690c0"
immune_col[3]<-"#4292c6"
immune_col[4]<-"#9ecae1"
immune_col[9]<-"#6baed6"
immune_col[21]<-"#9ecae1"
immune_col[22]<-"#6baed6"
show_col(immune_col)

cell_cols<-c(epi_col, immune_col, stromal_col)
names(cell_cols)<-c(epithelial, immune, stromal_and_glial)

ggplot(plt, aes(UMAP_1, UMAP_2, color=cluster_ID)) + geom_point(size=0.25) +
  theme_classic() +th_present +
  scale_color_manual(name="Cell",values = cell_cols, drop = T, guide=F)+
  xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-2, y=-4, label="Immune", color="#1C6BB0")+
  annotate("text",x=10, y=-10, label="Epithelial", color="#BB1419")+
  annotate("text",x=5, y=14, label="Stromal and Glial", color="#41AB5D")

ggsave(file=here("figs","primary_controls_UMAP_celltype_compartment.pdf"), w=5.5,h=5)
ggsave(file=here("figs/jpeg","primary_controls_UMAP_celltype_compartment.jpeg"), w=5.5,h=5)

####################################################################################################################

##############
## Expression all compartments
##############

MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")
crypt_villis = c("SEPP1", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8", "APOBEC1", "APOB", "APOA4", "APOA1", "NPC1L1", "EGFR", "KLF4", "ENPP3", "NT5E", "SLC28A2", "ADA")

## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))

## add cell type labels from split analysis
load(here("output","cell_label_whole_clustering.RData"))
identical(rownames(cell_label),rownames(d10x.primary@meta.data))
d10x.primary <- AddMetaData(d10x.primary, metadata = cell_label)

head(d10x.primary@meta.data)

d10x.primary.control<-subset(d10x.primary, subset = orig.ident == "Ctrl")
d10x.primary.control

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.primary.control <- NormalizeData(d10x.primary.control,scale.factor = 10000, normalization.method = "LogNormalize")
d10x.primary.control_scores<-FetchData(object = d10x.primary.control, vars = c(MHCI,crypt_villis))


### save MHC I score too
d10x.primary.control <- AddModuleScore(
  object = d10x.primary.control,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score'
)

d10x.primary.control <- AddModuleScore(
  object = d10x.primary.control,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)

score_data<-d10x.primary.control@meta.data[,c("MHCI_score1","crypt_villis_score1")]

d10x.exp.GOI<-FetchData(object = d10x.primary.control, vars = "NLRC5")


save(plt,d10x.primary.control_scores, score_data, d10x.exp.GOI, file=here("data","Primary_UMAP_controls_plus_MHCI_raw.RData"))


############
## UMAP Score
############
load(here("data","Primary_UMAP_controls_plus_MHCI_raw.RData"))
score_data$cell<-rownames(score_data)
plt_MHCIscore<-merge(plt,score_data, by="cell")

colours= brewer.pal(name="Blues", n=nlevels(diamonds$cut))
names(colours)= rev(levels(diamonds$cut))

ggplot(plt_MHCIscore, aes(UMAP_1, UMAP_2, color=MHCI_score1)) + geom_point(size=0.25) +
  theme_classic() +th_present +
  xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-2, y=-4, label="Immune", color="#1C6BB0")+
  annotate("text",x=10, y=-10, label="Epithelial", color="#BB1419")+
  annotate("text",x=5, y=14, label="Stromal and Glial", color="#41AB5D")+
  scale_color_viridis_c()

ggsave(file=here("figs","primary_controls_UMAP_celltype_compartment_MHCI.pdf"), w=5.5,h=5)
ggsave(file=here("figs/jpeg","primary_controls_UMAP_celltype_compartment_MHCI.jpeg"), w=5.5,h=5)

d10x.exp.GOI$cell<-rownames(d10x.exp.GOI)
plt<-merge(plt, d10x.exp.GOI, by="cell")
plt<-plt[(order(plt$NLRC5)),]

plot_grid(ggplot(plt, aes(UMAP_1, UMAP_2, color=cluster_ID)) + geom_point(size=0.25) +
            theme_classic() +th_present +
            scale_color_manual(name="Cell",values = cell_cols, drop = T, guide=F)+
            xlab("UMAP 1")+ylab("UMAP 2")+
            annotate("text",x=-2, y=-4, label="Immune", color="#1C6BB0")+
            annotate("text",x=10, y=-10, label="Epithelial", color="#BB1419")+
            annotate("text",x=7, y=14, label="Stromal and Glial", color="#41AB5D")+
            annotate("text",x=-11, y=-11, label="n = 27,622", color="grey40"),
          ggplot(plt_MHCIscore, aes(UMAP_1, UMAP_2, color=MHCI_score1)) + geom_point(size=0.25) +
            theme_classic() +th_present +
            xlab("UMAP 1")+ylab("UMAP 2")+
            annotate("text",x=-2, y=-4, label="Immune", color="#1C6BB0")+
            annotate("text",x=10, y=-10, label="Epithelial", color="#BB1419")+
            annotate("text",x=7, y=14, label="Stromal and Glial", color="#41AB5D")+
            annotate("text",x=-11, y=-11, label="n = 27,622", color="grey40")+
            scale_color_viridis_c(name="MHC I\nScore") + theme(legend.position = c(0.08, 0.8),
                                                               legend.text=element_text(size=8),
                                                               legend.title=element_text(size=10),
                                                               legend.key.height = unit(0.4, 'cm'),
                                                               legend.key.width = unit(0.2, 'cm')),
          ggplot(plt, aes(UMAP_1, UMAP_2, color=NLRC5))+geom_point(size=0.25)+
            th_present+theme_classic()+
            xlab("UMAP 1")+ylab("UMAP 2")+
            annotate("text",x=-2, y=-4, label="Immune", color="#1C6BB0")+
            annotate("text",x=10, y=-10, label="Epithelial", color="#BB1419")+
            annotate("text",x=7, y=14, label="Stromal and Glial", color="#41AB5D")+
            annotate("text",x=-11, y=-11, label="n = 27,622", color="grey40")+
            scale_color_continuous_sequential("ag_GrnYl",begin=0.0001, rev = F,na.value = "grey80")+theme(legend.position = c(0.08, 0.8),
                                                 legend.text=element_text(size=8),
                                                 legend.title=element_text(size=10),
                                                 legend.key.height = unit(0.4, 'cm'),
                                                 legend.key.width = unit(0.2, 'cm')),
          align="v",axis="lr", ncol=1)

ggsave(file=here("figs","primary_controls_UMAP_celltype_compartment_MHCI_both.pdf"), w=5.5,h=10.5)
ggsave(file=here("figs/jpeg","primary_controls_UMAP_celltype_compartment_MHCI_both.jpeg"), w=4.5,h=8.5)


# ggsave(file=here("../EBI/scRNAseq_codon/figs","primary_controls_UMAP_celltype_compartment_MHCI_both.pdf"), w=5.5,h=12.5)
# ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","primary_controls_UMAP_celltype_compartment_MHCI_both.jpeg"), w=5.5,h=12.5)




####################################################################################################################
## Correlation MHC I genes Controls only
####################################################################################################################
# 
# 
# MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")
# crypt_villis = c("SEPP1", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8", "APOBEC1", "APOB", "APOA4", "APOA1", "NPC1L1", "EGFR", "KLF4", "ENPP3", "NT5E", "SLC28A2", "ADA")
# 
# 
# ## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# # but not normalized or scaled
# d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))
# 
# ## add cell type labels from split analysis
# load(here("output","cell_label_whole_clustering.RData"))
# identical(rownames(cell_label),rownames(d10x.primary@meta.data))
# d10x.primary <- AddMetaData(d10x.primary, metadata = cell_label)
# 
# d10x.primary@meta.data$general_type<-"Immune"
# d10x.primary@meta.data$general_type[which(d10x.primary@meta.data$seurat_clusters%in%c(12,14,25))]<-"Epithelial"
# d10x.primary@meta.data$general_type[which(d10x.primary@meta.data$seurat_clusters%in%c(5,10,13,15,16,22,23,28,24))]<-"Stromal"
# 
# head(d10x.primary@meta.data)
# 
# d10x.primary_epi_raw<-subset(d10x.primary, subset = general_type == "Epithelial")
# d10x.primary_epi_raw<-subset(d10x.primary_epi_raw, subset = orig.ident == "Ctrl")
# d10x.primary_epi_raw
# 
# load(here("data","primary.epi.cells.RData"))
# cell_labels<-primary.epi.cells@meta.data
# 
# d10x.primary_epi_raw <- AddMetaData(d10x.primary_epi_raw, metadata = cell_labels)
# #filter remaining immune (179 cells)
# d10x.primary_epi_raw<-subset(d10x.primary_epi_raw, subset = cluster_ID != "B cell" & cluster_ID != "T cell")
# 
# 
# ##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# # This is log(TP10K+1)
# d10x.primary_epi_raw <- NormalizeData(d10x.primary_epi_raw,scale.factor = 10000, normalization.method = "LogNormalize")
# d10x.primary_epi_raw_scores<-FetchData(object = d10x.primary_epi_raw, vars = c(MHCI,crypt_villis))
# 


# 
# 
# ## diagional
# cormat<-cor(d10x.primary_epi_raw_scores, method="spearman")
# 
# fill_mat<-cormat
# fill_mat[upper.tri(fill_mat)]<- NA
# fill_melted_cormat <- melt(fill_mat, na.rm = TRUE)
# 
# ## remove same gene cor
# fill_melted_cormat<-fill_melted_cormat[which(!(fill_melted_cormat$Var1==fill_melted_cormat$Var2)),]
# 
# ## select some correlations to label
# label_max_offset<-fill_melted_cormat[which(fill_melted_cormat$Var1%in%crypt_villis & fill_melted_cormat$Var2%in%MHCI),]
# label_max_offset<-label_max_offset[!duplicated(label_max_offset),]
# 
# label_NLRC5<-fill_melted_cormat[which(fill_melted_cormat$Var1=="NLRC5" & fill_melted_cormat$Var2%in%MHCI),]
# 
# selected_labels<-rbind(label_NLRC5,label_max_offset[rev(order(label_max_offset$value)),][1:10,])
# 
# fill_melted_cormat$Var1<-factor(fill_melted_cormat$Var1, levels=c(MHCI,crypt_villis))
# fill_melted_cormat$Var2<-factor(fill_melted_cormat$Var2, levels=c(MHCI,crypt_villis))
# 
# 
# ggplot()+geom_tile(aes(Var1, Var2, fill=value),fill_melted_cormat)+
#   geom_text(aes(Var1, Var2, label=round(value,2)), selected_labels, color="grey40", size=2.5)+
#   scale_fill_gradient2(high="#2171b5",mid="#f3f7fc",low="red",
#                        na.value="grey90", midpoint=0)+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_text(size =9, color="black"),
#         legend.text = element_text(size =10),
#         legend.title = element_text(size =11),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())+
#   scale_y_discrete(position = "right") +
#   geom_rect(aes(xmin=length(MHCI)-0.5, xmax=length(c(MHCI,crypt_villis))-1.5,ymax=length(MHCI)+0.5, ymin=0.5),fill=NA, color="grey")
# 
# ggsave(file=here("figs/jpeg","score_correlation_primary_controls.jpeg"), w=10, h=9)
# ggsave(file=here("figs","score_correlation_primary_controls.pdf"), w=10, h=9)
# 
# 
# 
# #################
# ## Are these correlations more than expected by chance?
# #################
# label_max_offset<-fill_melted_cormat[which(fill_melted_cormat$Var1%in%crypt_villis & fill_melted_cormat$Var2%in%MHCI),]
# median(abs(label_max_offset$value))
# 
# 
# 
# #compare MHCI to random
# random_medians_MHCtornd<-sapply(1:1000, function(x){
#   set.seed(x)
#   ## random gene list
#   rnd<-sample(rownames(d10x.primary_epi_raw),length(crypt_villis))
#   
#   ## fetch random gene list
#   d10x.stim_MHCI_rnd<-FetchData(object = d10x.primary_epi_raw, vars = c(MHCI,rnd))
#   
#   ## diagional
#   cormat<-cor(d10x.stim_MHCI_rnd, method="spearman")
#   
#   fill_mat_rnd<-cormat
#   fill_mat_rnd[upper.tri(fill_mat_rnd)]<- NA
#   fill_melted_cormat_rnd <- melt(fill_mat_rnd, na.rm = TRUE)
#   
#   ## remove same gene cor
#   fill_melted_cormat_rnd<-fill_melted_cormat_rnd[which(!(fill_melted_cormat_rnd$Var1==fill_melted_cormat_rnd$Var2)),]
#   
#   ## select some correlations to label
#   label_max_offset_rnd<-fill_melted_cormat_rnd[which(fill_melted_cormat_rnd$Var1%in%rnd & fill_melted_cormat_rnd$Var2%in%MHCI),]
#   median(abs(label_max_offset_rnd$value))
# })
# 
# 
# ## how many larger
# length(which(random_medians_MHCtornd>=median(abs(label_max_offset$value))))
# ggplot()+geom_histogram(aes(random_medians_MHCtornd))+geom_vline(aes(xintercept=median(abs(label_max_offset$value))))
# 
# (length(which(random_medians_MHCtornd>=median(abs(label_max_offset$value))))+1)/(1000+1)
# 
# 
# 
# 
# 
# 
# #compare crypt villus to random
# random_medians_crypviltornd<-sapply(1:100, function(x){
#   set.seed(x)
#   ## random gene list
#   rnd<-sample(rownames(d10x.primary_epi_raw),length(MHCI))
#   
#   ## fetch random gene list
#   d10x.stim_crypt_villis_rnd<-FetchData(object = d10x.primary_epi_raw, vars = c(rnd,crypt_villis))
#   
#   ## diagional
#   cormat<-cor(d10x.stim_crypt_villis_rnd, method="spearman")
#   
#   fill_mat_rnd<-cormat
#   fill_mat_rnd[upper.tri(fill_mat_rnd)]<- NA
#   fill_melted_cormat_rnd <- melt(fill_mat_rnd, na.rm = TRUE)
#   
#   ## remove same gene cor
#   fill_melted_cormat_rnd<-fill_melted_cormat_rnd[which(!(fill_melted_cormat_rnd$Var1==fill_melted_cormat_rnd$Var2)),]
#   
#   ## select some correlations to label
#   label_max_offset_rnd<-fill_melted_cormat_rnd[which(fill_melted_cormat_rnd$Var1%in%crypt_villis & fill_melted_cormat_rnd$Var2%in%rnd),]
#   median(abs(label_max_offset_rnd$value))
# })
# 
# 
# ## how many larger
# length(which(random_medians_crypviltornd>=median(abs(label_max_offset$value))))
# ggplot()+geom_histogram(aes(random_medians_crypviltornd))+geom_vline(aes(xintercept=median(abs(label_max_offset$value))))
# 
# (length(which(random_medians_crypviltornd>=median(abs(label_max_offset$value))))+1)/(100+1)
# 
# 
# 
# 
# 
# ## random gene list
# rnd<-sample(rownames(d10x.primary_epi_raw),length(crypt_villis))
# 
# ## fetch random gene list
# d10x.stim_MHCI_rnd<-FetchData(object = d10x.primary_epi_raw, vars = c(MHCI,rnd))
# 
# ## diagional
# cormat<-cor(d10x.stim_MHCI_rnd)
# 
# fill_mat_rnd<-cormat
# fill_mat_rnd[upper.tri(fill_mat_rnd)]<- NA
# fill_melted_cormat_rnd <- melt(fill_mat_rnd, na.rm = TRUE)
# 
# ## remove same gene cor
# fill_melted_cormat_rnd<-fill_melted_cormat_rnd[which(!(fill_melted_cormat_rnd$Var1==fill_melted_cormat_rnd$Var2)),]
# 
# ## select some correlations to label
# label_max_offset_rnd<-fill_melted_cormat_rnd[which(fill_melted_cormat_rnd$Var1%in%rnd & fill_melted_cormat_rnd$Var2%in%MHCI),]
# median(abs(label_max_offset_rnd$value))
# 
# 
# label_NLRC5<-fill_melted_cormat_rnd[which(fill_melted_cormat_rnd$Var1=="NLRC5" & fill_melted_cormat_rnd$Var2%in%MHCI),]
# 
# selected_labels<-rbind(label_NLRC5,label_max_offset_rnd[rev(order(label_max_offset_rnd$value)),][1:10,])
# 
# fill_melted_cormat_rnd$Var1<-factor(fill_melted_cormat_rnd$Var1, levels=c(MHCI,rnd))
# fill_melted_cormat_rnd$Var2<-factor(fill_melted_cormat_rnd$Var2, levels=c(MHCI,rnd))
# 
# 
# ggplot()+geom_tile(aes(Var1, Var2, fill=value),fill_melted_cormat_rnd)+
#   geom_text(aes(Var1, Var2, label=round(value,2)), selected_labels, color="grey40", size=2.5)+
#   scale_fill_gradient2(high="#2171b5",mid="#f3f7fc",low="red",
#                        na.value="grey90", midpoint=0)+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_text(size =9, color="black"),
#         legend.text = element_text(size =10),
#         legend.title = element_text(size =11),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())+
#   scale_y_discrete(position = "right") +
#   geom_rect(aes(xmin=length(MHCI)-0.5, xmax=length(c(MHCI,rnd))-0.5,ymax=length(MHCI)+0.5, ymin=0.5),fill=NA, color="grey")
# 
# 
# #######################                  
# ## Individual correlation plots
# #######################
# 
# 
# gene_correlation<-function(gene_name1, gene_name2){
#   plt_gene<-FetchData(object = d10x.primary_epi_raw, vars = c(gene_name1, gene_name2))
#   #plt_gene$cell<-rownames(plt_gene)
#   
#   ## ignore 0s for correlation?
#   plt_gene$zero<-"Both non-zero"
#   plt_gene$zero[which(plt_gene[,1]==0 | plt_gene[,2]==0)]<-"One or both\ngenes zero"
#   
#   allcor<-cor(plt_gene[,1],plt_gene[,2], method="spearman")
#   #non-zero cor
#   plt_gene_non_zero<-plt_gene[which(plt_gene$zero=="Both non-zero"),]
#   nonzerocor<-cor(plt_gene_non_zero[,1],plt_gene_non_zero[,2], method="spearman")
#   
#   scatter<-ggplot()+geom_point(aes(plt_gene[,1], plt_gene[,2],  color=zero),plt_gene)+ylab(gene_name2)+xlab(gene_name1)+
#     theme_bw()+th+scale_color_manual(values=c("cornflowerblue","lightgrey"), name="Both genes\ndetected in cell")+
#     stat_smooth(aes(plt_gene[,1], plt_gene[,2]),plt_gene ,method="lm", se=F, color="black")+#
#     stat_smooth(aes(plt_gene_non_zero[,1], plt_gene_non_zero[,2]),plt_gene_non_zero ,method="lm", se=F, color="#1f66e5")+
#     theme(legend.text=element_text(size=8),
#           legend.title=element_text(size=10))
#   
#   get_leg = function(a.gplot){
#     tmp <- ggplot_gtable(ggplot_build(a.gplot))
#     leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#     legend <- tmp$grobs[[leg]]
#     legend}
#   
#   leg_scatter = get_leg(scatter)
#   
#   
#   density1<-ggplot(plt_gene, aes(plt_gene[,1]))+geom_density(fill="#a6bddb")+
#     theme_void()+th+xlab(gene_name1)+theme(axis.title=element_blank(),
#                                            axis.text=element_blank(),
#                                            axis.ticks=element_blank(),
#                                            legend.position = "none",
#                                            panel.background = element_blank(),
#                                            panel.grid.major = element_blank(), 
#                                            panel.grid.minor = element_blank(),
#                                            plot.margin = unit(c(0.75,0.2,-0.1,1), "cm"))
#   
#   density2<-ggplot(plt_gene, aes(plt_gene[,2]))+geom_density(fill="#a6bddb")+
#     coord_flip()+theme_void()+th+xlab(gene_name1)+theme(axis.title=element_blank(),
#                                                         axis.text=element_blank(),
#                                                         axis.ticks=element_blank(),
#                                                         legend.position = "none",
#                                                         panel.background = element_blank(),
#                                                         panel.grid.major = element_blank(), 
#                                                         panel.grid.minor = element_blank(),
#                                                         plot.margin = unit(c(0.2,1,1.1,-0.1), "cm"))
#   
#   placeholder<-ggplot() + theme_void()
#   
#   corvalues<-ggplot(plt_gene,aes(plt_gene[,1], plt_gene[,2],  color=zero)) + theme_void() +
#     annotate("text", x=max(plt_gene[,1])*0.25, y=max(plt_gene[,2])*0.5, label = paste("All cells\nR =",round(allcor,2)))+
#     annotate("text", x=max(plt_gene[,1])*0.75, y=max(plt_gene[,2])*0.5, label = paste("Both genes non-zero\nR =",round(nonzerocor,2)), color="#1f66e5")+
#     xlim(0,max(plt_gene[,1]))+ylim(0,max(plt_gene[,2]))
#   
#   grid_fig<-plot_grid(density1, leg_scatter, 
#                       scatter+theme(legend.position = "none"), density2, 
#                       corvalues,placeholder,
#                       ncol = 2, axis="lr", 
#                       rel_heights = c(1,4,0.5),
#                       rel_widths = c(4,1))
#   
#   grid_fig
#   
#   ggsave(file=paste(here("figs/jpeg/"),gene_name1,"_",gene_name2,"_correlation_primary_controls.jpeg", sep=""), w=7, h=7)
#   ggsave(file=paste(here("figs/"),gene_name1,"_",gene_name2,"_correlation_primary_controls.pdf", sep=""), w=7, h=7)
# }
# 
# gene_correlation("NLRC5","TAP1")
# gene_correlation("NLRC5","PSMB9")
# gene_correlation("B2M","HLA-B")
# gene_correlation("NLRC5","IRF1")
# 
# 
# gene_correlation("PLAC8","B2M")
# gene_correlation("IFI27","HLA-A")
# 
# 
# 
# ######### 
# ## summary statistics on correlation
# ######### 
# 
# gene_combos<-combn(MHCI, 2)
# 
# gene_correlation_values<-do.call(rbind, lapply(1:ncol(gene_combos), function(x){
#   gene_name1=gene_combos[1,x]
#   gene_name2=gene_combos[2,x]
#   plt_gene<-FetchData(object = d10x.primary_epi_raw, vars = c(gene_name1, gene_name2))
#   
#   ## ignore 0s for correlation
#   plt_gene$zero<-"Both non-zero"
#   plt_gene$zero[which(plt_gene[,1]==0 | plt_gene[,2]==0)]<-"One or both\ngenes zero"
#   
#   allcor<-cor(plt_gene[,1],plt_gene[,2], method="spearman")
#   #non-zero cor
#   plt_gene_non_zero<-plt_gene[which(plt_gene$zero=="Both non-zero"),]
#   nonzerocor<-cor(plt_gene_non_zero[,1],plt_gene_non_zero[,2], method="spearman")
#   data.frame(gene_name1=gene_name1, gene_name2=gene_name2, nonzerocor=nonzerocor, allcor=allcor)
# }))
# 
# gene_correlation_values$nonzerocor[which(is.na(gene_correlation_values$nonzerocor))]<-0
# 
# 
# 
# 
# plt_gene_MHCI<-FetchData(object = d10x.primary_epi_raw, vars = MHCI)
# counts_cells<-apply(plt_gene_MHCI, 2, function(x) length(which(x!=0)))
# mean_cells<-apply(plt_gene_MHCI, 2, mean)
# 
# node_values<-cbind(melt(counts_cells), melt(mean_cells))
# colnames(node_values)<-c("cell_expressing","mean_expression")
# node_values$gene=rownames(node_values)
# 
# 
# 
# ##### MHCI correlation only in non zero alternate to network
# gene_correlation_values$gene_name1<-factor(gene_correlation_values$gene_name1, levels=rev(MHCI))
# gene_correlation_values$gene_name2<-factor(gene_correlation_values$gene_name2, levels=rev(MHCI))
# 
# selected_labels<-gene_correlation_values[which(gene_correlation_values$gene_name2=="NLRC5"),]
# 
# ggplot()+geom_tile(aes(gene_name1, gene_name2, fill=nonzerocor),gene_correlation_values, color="black")+
#   geom_text(aes(gene_name1, gene_name2, label=round(nonzerocor,2)), selected_labels, color="grey10", size=4)+
#   scale_fill_gradient2(high="#2171b5",mid="#f3f7fc",low="red",
#                        na.value="grey90", midpoint=0.1, name="Correlation", breaks=seq(-1,1,0.2), limits=c(-1,1))+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         axis.ticks = element_blank(),
#         axis.title = element_blank(),
#         axis.text = element_text(size =12, color="black"),
#         legend.text = element_text(size =10),
#         legend.title = element_text(size =11),
#         panel.background = element_blank(),
#         panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank())+
#   scale_y_discrete(position = "right") 
# 
# ggsave(file=here("figs/jpeg","MHCI_score_correlation_primary_nonzero_controls.jpeg"), w=10, h=9)
# ggsave(file=here("figs","MHCI_score_correlation_primary_nonzero_controls.pdf"), w=10, h=9)
# 


###individual gene expression
load(here("data/Primary_UMAP_controls_plus_MHCI_raw.RData"))

d10x.primary.control_scores$cell<-rownames(d10x.primary.control_scores)
plt$general_type<-as.factor(plt$general_type)
levels(plt$general_type)<-c("Epithelial","Immune","Stromal\nand\nGlial")
comp_simple<-list(c("Epithelial","Immune"),c("Epithelial","Stromal\nand\nGlial"),c("Immune","Stromal\nand\nGlial"))


plt_gene_compartment<-function(gene){
  d10x.primary.control_scores_gene<-d10x.primary.control_scores[, c("cell", gene)]
  colnames(d10x.primary.control_scores_gene)<-c("cell","gene_name")
  plt_gene<-merge(plt,d10x.primary.control_scores_gene, by="cell")
  
  print(pairwise.t.test(plt_gene$gene_name, plt_gene$general_type,p.adjust.method = "BH", pool.sd = FALSE))
  ggplot(plt_gene, aes(general_type,gene_name))+
    geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=general_type))+
    theme_bw()+th_present+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
    scale_fill_manual(values=c("#b2182b","#4393c3","#5aae61"), guide=F)+ylab(paste(gene, "Expression"))+xlab("Compartment")+
    geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
                size = 0.3,vjust = 0.5,
                textsize = 3,  map_signif_level = T, color="grey60")#
}

plt_gene_compartment("B2M")#
plt_gene_compartment("IRF1")

plt_gene_compartment("NLRC5")
ggsave(file=here("figs","primary_controls_compartment_NLRC5.pdf"), w=5,h=7)
ggsave(file=here("figs/jpeg","primary_controls_compartment_NLRC5.jpeg"), w=3,h=6)


plt_gene_compartment("TAP1")#
plt_gene_compartment("TAP2")#

plt_gene_compartment("HLA-A")
plt_gene_compartment("HLA-B")
plt_gene_compartment("HLA-C")#
plt_gene_compartment("HLA-E")#
plt_gene_compartment("HLA-F")#
plt_gene_compartment("HLA-G")

plt_gene_compartment("PSMB9")
plt_gene_compartment("PSMB8")#


MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")

plot_grid(plotlist=lapply(MHCI,plt_gene_compartment))
ggsave(file=here("figs","primary_controls_compartment_MHCI.pdf"), w=15,h=20)
ggsave(file=here("figs/jpeg","primary_controls_compartment_MHCI.jpeg"), w=15,h=20)



plt_gene_compartment_log<-function(gene){
  d10x.primary.control_scores_gene<-d10x.primary.control_scores[, c("cell", gene)]
  colnames(d10x.primary.control_scores_gene)<-c("cell","gene_name")
  plt_gene<-merge(plt,d10x.primary.control_scores_gene, by="cell")
  
  print(pairwise.t.test(plt_gene$gene_name, plt_gene$general_type,p.adjust.method = "BH", pool.sd = FALSE))
  ggplot(plt_gene, aes(general_type,log(gene_name)))+
    geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=general_type))+
    theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
    scale_fill_manual(values=c("#b2182b","#4393c3","#5aae61"), guide=F)+ylab(paste(gene, "Expression"))+xlab("Compartment")+
    geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
                size = 0.3,vjust = 0.5,
                textsize = 3,  map_signif_level = T, color="grey60")#
}


plt_gene_compartment_log("B2M")
plt_gene_compartment_log("IRF1")

plt_gene_compartment_log("NLRC5")
plt_gene_compartment_log("TAP1")
plt_gene_compartment_log("TAP2")

plt_gene_compartment_log("HLA-A")
plt_gene_compartment_log("HLA-B")
plt_gene_compartment_log("HLA-C")
plt_gene_compartment_log("HLA-E")
plt_gene_compartment_log("HLA-F")
plt_gene_compartment_log("HLA-G")

plt_gene_compartment_log("PSMB9")
plt_gene_compartment_log("PSMB8")


###########
## control only MHCI
###########
score_data$cell<-rownames(score_data)
plt_MHCIscore<-merge(plt,score_data, by="cell")

print(pairwise.t.test(plt_MHCIscore$MHCI_score1, plt_MHCIscore$general_type,p.adjust.method = "BH", pool.sd = FALSE))
print(pairwise.t.test(plt$NLRC5, plt$general_type,p.adjust.method = "BH", pool.sd = FALSE))

plot_grid(
ggplot(plt_MHCIscore, aes(general_type,MHCI_score1))+
  geom_violin(fill="grey80", color="grey80")+geom_boxplot(width=0.1,aes(fill=general_type),outlier.shape = NA)+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12))+
  scale_fill_manual(values=c("#b2182b","#4393c3","#5aae61"), guide=F)+ylab("MHC I Score")+xlab("")+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 2,  map_signif_level = T, color="grey60"),
ggplot(plt, aes(general_type,NLRC5))+
  geom_violin(fill="grey80", color="grey80")+geom_boxplot(width=0.1,aes(fill=general_type),outlier.shape = NA)+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12))+
  scale_fill_manual(values=c("#b2182b","#4393c3","#5aae61"), guide=F)+ylab("NLRC5 Count")+xlab("")+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 2,  map_signif_level = T, color="grey60"))

ggsave(file=here("figs","primary_controls_boxplot_compartment_MHCI_both.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","primary_controls_boxplot_compartment_MHCI_both.jpeg"), w=5,h=4)

# ggsave(file=here("../EBI/scRNAseq_codon/figs","primary_controls_boxplot_compartment_MHCI_both.pdf"), w=5,h=4)
# ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","primary_controls_boxplot_compartment_MHCI_both.jpeg"), w=5,h=4)
# 
