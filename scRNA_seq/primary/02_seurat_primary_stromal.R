#'---
#'title: scRNAseq primary epithelial clustering
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


## load raw primay data 
d10x.primary.raw<-readRDS(here("data","d10x_primary_raw_merged.rds"))
d10x.primary.raw

## pull comparment assignments from combined primary clustering/marker identification
load(here("output","cell_label_whole_clustering.RData"))


## add cell label information
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x.primary.raw), cell_label$index),]
identical(colnames(d10x.primary.raw), cell_label$index)

d10x.primary.raw<- AddMetaData(d10x.primary.raw, cell_label$seurat_clusters, col.name = "seurat_wholePrimary_clusters")
d10x.primary.raw<- AddMetaData(d10x.primary.raw, cell_label$elmentaite_annotation_V2, col.name = "elmentaite_annotation_V2")

d10x.primary.raw@meta.data$general_type<-"Immune"
d10x.primary.raw@meta.data$general_type[which(d10x.primary.raw@meta.data$seurat_wholePrimary_clusters%in%c(12,14,25))]<-"Epithelial"
d10x.primary.raw@meta.data$general_type[which(d10x.primary.raw@meta.data$seurat_wholePrimary_clusters%in%c(5,10,13,15,16,22,23,28,24))]<-"Stromal"

table(d10x.primary.raw@meta.data$general_type)

## subset to just stromal
primary.stromal.cells<-subset(d10x.primary.raw, subset = general_type == "Stromal")

## remove neonatal
primary.stromal.cells<-subset(primary.stromal.cells, subset = orig.ident != "neo")

# An object of class Seurat
# 24945 features across 45612 samples within 1 assay
# Active assay: RNA (24945 features, 0 variable features)


primary.stromal.cells <- NormalizeData(primary.stromal.cells)
primary.stromal.cells <- FindVariableFeatures(primary.stromal.cells, selection.method = "vst", nfeatures = 2000)
primary.stromal.cells <- ScaleData(primary.stromal.cells) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
primary.stromal.cells <- RunPCA(primary.stromal.cells, ndims.print = 1:10, nfeatures.print = 10)
primary.stromal.cells <- RunUMAP(primary.stromal.cells, dims = 1:30)
primary.stromal.cells <- RunTSNE(primary.stromal.cells, dims = 1:30)



######################
## cell cycle gene expression
######################

## which PC associated to nfeature? 2
pca_mat<-as.data.frame(Embeddings(object = primary.stromal.cells, reduction = "pca"))
pca_mat$cell<-rownames(pca_mat)
meta_primary<-primary.stromal.cells@meta.data
meta_primary$cell<-rownames(meta_primary)
plt<-merge(pca_mat, meta_primary,by="cell")

sapply(2:51, function(x) cor(plt[,x],plt[,"nFeature_RNA"]))
FeaturePlot(primary.stromal.cells, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
## PC2 tops genes include cell cycle genes: "TPX2"   "TOP2A"  "MKI67"  "CENPF"  "SMC4"   "TUBB4B"

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

primary.stromal.cells <- CellCycleScoring(primary.stromal.cells, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(primary.stromal.cells, reduction="pca")
ggsave(file=here("figs","stromal_primary_cell_cycle_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","stromal_primary_cell_cycle_PCA.jpeg"), w=5,h=5)

FeaturePlot(primary.stromal.cells, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","stromal_primary_nfeature_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","stromal_primary_nfeature_PCA.jpeg"), w=5.25,h=5)

DimPlot(primary.stromal.cells, reduction="umap", group.by = "Phase")


## regress out cell cycle and other covariates
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#By default, sctransform accounts for cellular sequencing depth, or nUMIs.
primary.stromal.cells <- SCTransform(primary.stromal.cells, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE)


# dimension reduction
primary.stromal.cells <- RunPCA(primary.stromal.cells, verbose = FALSE)
primary.stromal.cells <- RunUMAP(primary.stromal.cells, dims = 1:30)
primary.stromal.cells <- RunTSNE(primary.stromal.cells, dims = 1:30)

# cluster
primary.stromal.cells <- FindNeighbors(primary.stromal.cells, reduction = "pca", dims = 1:20)
primary.stromal.cells <- FindClusters(primary.stromal.cells, resolution = 0.8)

# relabel meta data
primary.stromal.cells@meta.data$orig.ident[which(primary.stromal.cells@meta.data$orig.ident=="Ctrl")]<-"Control"
primary.stromal.cells@meta.data$orig.ident[which(primary.stromal.cells@meta.data$orig.ident=="neo")]<-"Neonatal"

saveRDS(primary.stromal.cells, file = here("data","d10x_stromal_primary_normalized.rds"))
print(table(primary.stromal.cells@meta.data$elmentaite_annotation_V2))







primary.stromal.cells<-readRDS(here("data","d10x_stromal_primary_normalized.rds"))


## visualize
DimPlot(primary.stromal.cells, reduction = "umap", group.by = "Phase", pt.size=1)
DimPlot(primary.stromal.cells, reduction = "umap", group.by = "Phase",split.by = "orig.ident", pt.size=1)


DimPlot(primary.stromal.cells, reduction = "umap", split.by = "orig.ident",label=T, pt.size=0.25)
ggsave(file=here("figs","stromal_primary_UMAP_source.pdf"), w=20,h=4)
ggsave(file=here("figs/jpeg","stromal_primary_UMAP_source.jpeg"), w=10,h=4)

DimPlot(primary.stromal.cells, reduction = "tsne", split.by = "orig.ident",label=T, pt.size=0.25)
ggsave(file=here("figs","stromal_primary_TSNE_raw_source.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","stromal_primary_TSNE__source.jpeg"), w=7,h=6)

DimPlot(primary.stromal.cells, reduction = "tsne",  pt.size=0.25)
ggsave(file=here("figs","stromal_primary_TSNE_raw_cluster.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","stromal_primary_TSNE_raw_cluster.jpeg"), w=7,h=6)

DimPlot(primary.stromal.cells, reduction = "umap",  pt.size=0.25, label=T)
ggsave(file=here("figs","stromal_primary_UMAP_raw_cluster.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","stromal_primary_UMAPraw_cluster.jpeg"), w=7,h=6)


DimPlot(primary.stromal.cells, reduction="pca", pt.size=1)

DimPlot(primary.stromal.cells, reduction = "umap", group.by = "individual", pt.size=1)
DimPlot(primary.stromal.cells, reduction = "umap", split.by = "individual", group.by = "orig.ident", pt.size=1)




table(primary.stromal.cells@meta.data$elmentaite_annotation_V2)

primary.stromal.cells@meta.data$elmentaite_annotation_V2<-as.factor(primary.stromal.cells@meta.data$elmentaite_annotation_V2)
primary.stromal.cells@meta.data$elmentaite_annotation_V2<-factor(primary.stromal.cells@meta.data$elmentaite_annotation_V2, levels = c(
  "Activated B cell","activated DC" ,"Activated T","Arterial endothelial cell" ,"B cell","CD4 T cell",
  "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
  "IgA plasma cell","IgG plasma cell","Lymphatic endothelial cell","Macrophage" ,"mast cells","Memory B cell" ,
  "Monocyte","pDC" ,"pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"TA","Tfh" ,
  "Treg", "Venous endothelial cell","Glial cell","myofibroblast","unknown_filtered"))

#' ### dev cell labels
cols_manual<-c(  "#238443","#810f7c" ,"#02818a","#8c510a" ,"#78c679","#67a9cf",
                 "#3690c0","#8073ac","#b2abd2","#238443","#dd3497","#1d91c0","#addd8e","#014636",
                 "#253494","#081d58","#bf812d","#7a0177" ,"#ce1256","#006d2c" ,
                 "#a50f15","#542788" ,"#e31a1c","#e08214","#ef6548","#fdd49e" ,"#f46d43","#4575b4" ,
                 "#016c59", "#bf812d","#c51b7d","#fb9a99","lightgrey")

names(cols_manual) <- levels(primary.stromal.cells@meta.data$elmentaite_annotation_V2)
fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)





DimPlot(primary.stromal.cells, reduction = "umap", group.by = "elmentaite_annotation_V2", pt.size=1, label=T)+colcale_cols_manual
ggsave(file=here("figs","primary_stromal_UMAP_cellLabel_devcell.pdf"), w=10,h=5)
ggsave(file=here("figs/jpeg","primary_stromal_UMAP_cellLabel_devcell.jpeg"), w=10,h=5)

DimPlot(primary.stromal.cells, reduction = "tsne", group.by = "elmentaite_annotation_V2", pt.size=1, label=T)+colcale_cols_manual
ggsave(file=here("figs","primary_stromal_tsne_cellLabel_devcell.pdf"), w=10,h=5)
ggsave(file=here("figs/jpeg","primary_stromal_tsne_cellLabel_devcell.jpeg"), w=10,h=5)

## plot know cells on top of unknown (UMAP)
umap_mat<-as.data.frame(Embeddings(object = primary.stromal.cells, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta_primary<-primary.stromal.cells@meta.data
meta_primary$cell<-rownames(meta_primary)

plt<-merge(umap_mat, meta_primary, by="cell")

primary.stromal.cells_known<-plt[which(plt$elmentaite_annotation_V2 != "unknown_filtered"),]
primary.stromal.cells_unknown<-plt[which(plt$elmentaite_annotation_V2 == "unknown_filtered"),]

ggplot()+
  geom_point(aes(UMAP_1,UMAP_2, color=elmentaite_annotation_V2),primary.stromal.cells_unknown,size=1.5)+
  geom_point(aes(UMAP_1,UMAP_2, color=elmentaite_annotation_V2),primary.stromal.cells_known,size=1.5)+
  colcale_cols_manual+theme_classic()+
  th+theme(plot.title = element_text(size = 16, face = "bold"))

ggsave(file=here("figs","primary_stromal_UMAP_cellLabel_devcell_nolab.pdf"), w=10,h=5)
ggsave(file=here("figs/jpeg","primary_stromal_UMAP_cellLabel_devcell_nolab.jpeg"), w=10,h=5)


## plot know cells on top of unknown (TSNE)
tsne_mat<-as.data.frame(Embeddings(object = primary.stromal.cells, reduction = "tsne"))#
tsne_mat$cell<-rownames(tsne_mat)

meta_primary<-primary.stromal.cells@meta.data
meta_primary$cell<-rownames(meta_primary)

plt<-merge(tsne_mat, meta_primary, by="cell")

primary.stromal.cells_known<-plt[which(plt$elmentaite_annotation_V2 != "unknown_filtered"),]
primary.stromal.cells_unknown<-plt[which(plt$elmentaite_annotation_V2 == "unknown_filtered"),]

ggplot()+
  geom_point(aes(tSNE_1,tSNE_2, color=elmentaite_annotation_V2),primary.stromal.cells_unknown,size=1.5)+
  geom_point(aes(tSNE_1,tSNE_2, color=elmentaite_annotation_V2),primary.stromal.cells_known,size=1.5)+
  colcale_cols_manual+theme_classic()+
  th+theme(plot.title = element_text(size = 16, face = "bold"))

ggsave(file=here("figs","primary_stromal_tsne_cellLabel_devcell_nolab.pdf"), w=10,h=5)
ggsave(file=here("figs/jpeg","primary_stromal_tsne_cellLabel_devcell_nolab.jpeg"), w=10,h=5)


DimPlot(primary.stromal.cells, reduction = "umap", pt.size=1, label=T)
ggsave(file=here("figs","primary_stromal_UMAP_cluster.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","primary_stromal_UMAP_cluster.jpeg"), w=7,h=5)

DimPlot(primary.stromal.cells, reduction = "tsne", pt.size=1, label=T)
ggsave(file=here("figs","primary_stromal_tsne_cluster.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","primary_stromal_tsne_cluster.jpeg"), w=7,h=5)



## relable clusters with cell types
table(primary.stromal.cells@meta.data$elmentaite_annotation_V2,primary.stromal.cells@meta.data$seurat_clusters)


primary.stromal.cells$cluster_ID<-"Unidentified"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==0)]<-"S1 fibroblasts"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==1)]<-"Arterial endothelial cell"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==2)]<-"S2 fibroblasts"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==3)]<-"S1 fibroblasts"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==4)]<-"S1 fibroblasts"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==5)]<-"pericyte"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==6)]<-"S4 fibroblasts"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==7)]<-"Venous endothelial cell"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==8)]<-"S4 fibroblasts"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==9)]<-"S4 fibroblasts"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==10)]<-"S4 fibroblasts"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==11)]<-"S1 fibroblasts"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==12)]<-"Lymphatic endothelial cell"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==13)]<-"Glial cell"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==14)]<-"Venous endothelial cell"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==15)]<-"S1 fibroblasts"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==16)]<-"pericyte"
primary.stromal.cells$cluster_ID[which(primary.stromal.cells$seurat_clusters==17)]<-"S2 fibroblasts"



DimPlot(primary.stromal.cells, reduction = "umap", group.by = "cluster_ID", pt.size=1, label=T)+colcale_cols_manual
ggsave(file=here("figs","primary_stromal_UMAP_label_expanded.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","primary_stromal_UMAP_label_expanded.jpeg"), w=7,h=5)


DimPlot(primary.stromal.cells, reduction = "umap", split.by = "orig.ident",  group.by = "cluster_ID",label=F, pt.size=0.25)+colcale_cols_manual
ggsave(file=here("figs","stromal_primary_UMAP_source_cluster.pdf"), w=10,h=4)
ggsave(file=here("figs/jpeg","stromal_primary_UMAP_source_cluster.jpeg"), w=10,h=4)


save(primary.stromal.cells, file=here("data","primary.stromal.cells.RData"))


stromal_cell_labels<-primary.stromal.cells@meta.data[,c("seurat_clusters","cluster_ID")]
save(stromal_cell_labels, file=here("output","strom_celltype_label.Rdata"))
