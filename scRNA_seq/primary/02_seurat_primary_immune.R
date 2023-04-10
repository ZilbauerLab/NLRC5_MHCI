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

## subset to just immune
primary.immune.cells<-subset(d10x.primary.raw, subset = general_type == "Immune")

## remove neonatal
primary.immune.cells<-subset(primary.immune.cells, subset = orig.ident != "neo")


# An object of class Seurat
# 24945 features across 45612 samples within 1 assay
# Active assay: RNA (24945 features, 0 variable features)


primary.immune.cells <- NormalizeData(primary.immune.cells)
primary.immune.cells <- FindVariableFeatures(primary.immune.cells, selection.method = "vst", nfeatures = 2000)
primary.immune.cells <- ScaleData(primary.immune.cells) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
primary.immune.cells <- RunPCA(primary.immune.cells, ndims.print = 1:10, nfeatures.print = 10)
primary.immune.cells <- RunUMAP(primary.immune.cells, dims = 1:30)
primary.immune.cells <- RunTSNE(primary.immune.cells, dims = 1:30)



######################
## cell cycle gene expression
######################

## which PC associated to nfeature? 2
pca_mat<-as.data.frame(Embeddings(object = primary.immune.cells, reduction = "pca"))
pca_mat$cell<-rownames(pca_mat)
meta_primary<-primary.immune.cells@meta.data
meta_primary$cell<-rownames(meta_primary)
plt<-merge(pca_mat, meta_primary,by="cell")

sapply(2:51, function(x) cor(plt[,x],plt[,"nFeature_RNA"]))
FeaturePlot(primary.immune.cells, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
## PC2 tops genes include cell cycle genes: "TPX2"   "TOP2A"  "MKI67"  "CENPF"  "SMC4"   "TUBB4B"

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

primary.immune.cells <- CellCycleScoring(primary.immune.cells, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(primary.immune.cells, reduction="pca")
ggsave(file=here("figs","immune_primary_cell_cycle_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","immune_primary_cell_cycle_PCA.jpeg"), w=5,h=5)

FeaturePlot(primary.immune.cells, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","immune_primary_nfeature_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","immune_primary_nfeature_PCA.jpeg"), w=5.25,h=5)

DimPlot(primary.immune.cells, reduction="umap", group.by = "Phase")


## regress out cell cycle and other covariates
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#By default, sctransform accounts for cellular sequencing depth, or nUMIs.
primary.immune.cells <- SCTransform(primary.immune.cells, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE)


# dimension reduction
primary.immune.cells <- RunPCA(primary.immune.cells, verbose = FALSE)
primary.immune.cells <- RunUMAP(primary.immune.cells, dims = 1:30)
primary.immune.cells <- RunTSNE(primary.immune.cells, dims = 1:30)

# cluster
primary.immune.cells <- FindNeighbors(primary.immune.cells, reduction = "pca", dims = 1:20)
primary.immune.cells <- FindClusters(primary.immune.cells, resolution = 0.8)

# relabel meta data
primary.immune.cells@meta.data$orig.ident[which(primary.immune.cells@meta.data$orig.ident=="Ctrl")]<-"Control"
primary.immune.cells@meta.data$orig.ident[which(primary.immune.cells@meta.data$orig.ident=="neo")]<-"Neonatal"

saveRDS(primary.immune.cells, file = here("data","d10x_immune_primary_normalized.rds"))






primary.immune.cells<-readRDS(here("data","d10x_immune_primary_normalized.rds"))



## visualize
DimPlot(primary.immune.cells, reduction = "umap", group.by = "Phase", pt.size=1)
DimPlot(primary.immune.cells, reduction = "umap", group.by = "Phase",split.by = "orig.ident", pt.size=1)


DimPlot(primary.immune.cells, reduction = "umap", split.by = "orig.ident",label=T, pt.size=0.25)
ggsave(file=here("figs","immune_primary_UMAP_source.pdf"), w=20,h=4)
ggsave(file=here("figs/jpeg","immune_primary_UMAP_source.jpeg"), w=10,h=4)

DimPlot(primary.immune.cells, reduction = "tsne", split.by = "orig.ident",label=T, pt.size=0.25)
ggsave(file=here("figs","immune_primary_TSNE_raw_source.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","immune_primary_TSNE__source.jpeg"), w=7,h=6)

DimPlot(primary.immune.cells, reduction = "tsne",  pt.size=0.25)
ggsave(file=here("figs","immune_primary_TSNE_raw_cluster.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","immune_primary_TSNE_raw_cluster.jpeg"), w=7,h=6)

DimPlot(primary.immune.cells, reduction = "umap",  pt.size=0.25, label=T)
ggsave(file=here("figs","immune_primary_UMAP_raw_cluster.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","immune_primary_UMAPraw_cluster.jpeg"), w=7,h=6)


DimPlot(primary.immune.cells, reduction="pca", pt.size=1)

DimPlot(primary.immune.cells, reduction = "umap", group.by = "individual", pt.size=1)
DimPlot(primary.immune.cells, reduction = "umap", split.by = "individual", group.by = "orig.ident", pt.size=1)




table(primary.immune.cells@meta.data$elmentaite_annotation_V2)

primary.immune.cells@meta.data$elmentaite_annotation_V2<-as.factor(primary.immune.cells@meta.data$elmentaite_annotation_V2)
primary.immune.cells@meta.data$elmentaite_annotation_V2<-factor(primary.immune.cells@meta.data$elmentaite_annotation_V2, levels = c(
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


names(cols_manual) <- levels(primary.immune.cells@meta.data$elmentaite_annotation_V2)
fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)





DimPlot(primary.immune.cells, reduction = "umap", group.by = "elmentaite_annotation_V2", pt.size=1, label=T)+colcale_cols_manual
ggsave(file=here("figs","primary_immune_UMAP_cellLabel_devcell.pdf"), w=10,h=5)
ggsave(file=here("figs/jpeg","primary_immune_UMAP_cellLabel_devcell.jpeg"), w=10,h=5)

DimPlot(primary.immune.cells, reduction = "tsne", group.by = "elmentaite_annotation_V2", pt.size=1, label=T)+colcale_cols_manual
ggsave(file=here("figs","primary_immune_tsne_cellLabel_devcell.pdf"), w=10,h=5)
ggsave(file=here("figs/jpeg","primary_immune_tsne_cellLabel_devcell.jpeg"), w=10,h=5)

## plot know cells on top of unknown (UMAP)
umap_mat<-as.data.frame(Embeddings(object = primary.immune.cells, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta_primary<-primary.immune.cells@meta.data
meta_primary$cell<-rownames(meta_primary)

plt<-merge(umap_mat, meta_primary, by="cell")

primary.immune.cells_known<-plt[which(plt$elmentaite_annotation_V2 != "unknown_filtered"),]
primary.immune.cells_unknown<-plt[which(plt$elmentaite_annotation_V2 == "unknown_filtered"),]

ggplot()+
  geom_point(aes(UMAP_1,UMAP_2, color=elmentaite_annotation_V2),primary.immune.cells_unknown,size=1.5)+
  geom_point(aes(UMAP_1,UMAP_2, color=elmentaite_annotation_V2),primary.immune.cells_known,size=1.5)+
  colcale_cols_manual+theme_classic()+
  th+theme(plot.title = element_text(size = 16, face = "bold"))

ggsave(file=here("figs","primary_immune_UMAP_cellLabel_devcell_nolab.pdf"), w=10,h=5)
ggsave(file=here("figs/jpeg","primary_immune_UMAP_cellLabel_devcell_nolab.jpeg"), w=10,h=5)


## plot know cells on top of unknown (TSNE)
tsne_mat<-as.data.frame(Embeddings(object = primary.immune.cells, reduction = "tsne"))#
tsne_mat$cell<-rownames(tsne_mat)

meta_primary<-primary.immune.cells@meta.data
meta_primary$cell<-rownames(meta_primary)

plt<-merge(tsne_mat, meta_primary, by="cell")

primary.immune.cells_known<-plt[which(plt$elmentaite_annotation_V2 != "unknown_filtered"),]
primary.immune.cells_unknown<-plt[which(plt$elmentaite_annotation_V2 == "unknown_filtered"),]

ggplot()+
  geom_point(aes(tSNE_1,tSNE_2, color=elmentaite_annotation_V2),primary.immune.cells_unknown,size=1.5)+
  geom_point(aes(tSNE_1,tSNE_2, color=elmentaite_annotation_V2),primary.immune.cells_known,size=1.5)+
  colcale_cols_manual+theme_classic()+
  th+theme(plot.title = element_text(size = 16, face = "bold"))

ggsave(file=here("figs","primary_immune_tsne_cellLabel_devcell_nolab.pdf"), w=10,h=5)
ggsave(file=here("figs/jpeg","primary_immune_tsne_cellLabel_devcell_nolab.jpeg"), w=10,h=5)


DimPlot(primary.immune.cells, reduction = "umap", pt.size=1, label=T)
ggsave(file=here("figs","primary_immune_UMAP_cluster.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","primary_immune_UMAP_cluster.jpeg"), w=7,h=5)

DimPlot(primary.immune.cells, reduction = "tsne", pt.size=1, label=T)
ggsave(file=here("figs","primary_immune_tsne_cluster.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","primary_immune_tsne_cluster.jpeg"), w=7,h=5)



## relable clusters with cell types
table(primary.immune.cells@meta.data$elmentaite_annotation_V2,primary.immune.cells@meta.data$seurat_clusters)


primary.immune.cells$cluster_ID<-"Unidentified"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==0)]<-"B cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==1)]<-"T cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==2)]<-"T cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==3)]<-"B cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==4)]<-"B cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==5)]<-"T cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==6)]<-"B cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==7)]<-"B cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==8)]<-"T cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==9)]<-"T cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==10)]<-"Plasma cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==11)]<-"Myeloid"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==12)]<-"T cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==13)]<-"Myeloid"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==14)]<-"T cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==15)]<-"Plasma cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==16)]<-"myofibroblast"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==17)]<-"Plasma cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==18)]<-"Plasma cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==19)]<-"Plasma cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==20)]<-"Mast cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==21)]<-"Plasma cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==22)]<-"Plasma cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==23)]<-"B cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==24)]<-"Plasma cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==25)]<-"S4 fibroblasts"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==26)]<-"Plasma cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==27)]<-"Plasma cell"
primary.immune.cells$cluster_ID[which(primary.immune.cells$seurat_clusters==28)]<-"Myeloid"


DimPlot(primary.immune.cells, reduction = "umap", pt.size=1, label=T)

DimPlot(primary.immune.cells, reduction = "umap", group.by="cluster_ID", pt.size=1, label=T)+scale_color_manual(values=c("#41ab5d","#fb6a4a","#ef3b2c","#4a1486","#4292c6","pink","blue"))
ggsave(file=here("figs","primary_immune_UMAP_cluster_broad_label.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","primary_immune_UMAP_cluster_broad_label.jpeg"), w=7,h=5)

primary.mast.cells<-subset(primary.immune.cells, subset = cluster_ID == "Mast cell")
primary.S4.cells<-subset(primary.immune.cells, subset = cluster_ID == "S4 fibroblasts")
primary.myofibroblast.cells<-subset(primary.immune.cells, subset = cluster_ID == "myofibroblast")


FeaturePlot(primary.immune.cells, features = "FCER2",reduction = "umap", min.cutoff = "q9", pt.size=1)
FeaturePlot(primary.immune.cells, features = "TNFRSF13B",reduction = "umap", min.cutoff = "q9", pt.size=1)
FeaturePlot(primary.immune.cells, features = "CD19",reduction = "umap", min.cutoff = "q9", pt.size=1)





#######################
### Next labelling iteration
######################
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
d10x.primary.raw@meta.data$general_type[which(d10x.primary.raw@meta.data$seurat_wholePrimary_clusters%in%c(5,10,13,15,16,22,23,28))]<-"Stromal"

table(d10x.primary.raw@meta.data$general_type)

## subset to broad catagories
primary.immune.cells.raw<-subset(d10x.primary.raw, subset = general_type == "Immune")
primary.immune.cells.raw<- AddMetaData(primary.immune.cells.raw, primary.immune.cells$cluster_ID, col.name = "broad_label")

primary.B.cells.raw<-subset(primary.immune.cells.raw, subset = broad_label == "B cell")
primary.T.cells.raw<-subset(primary.immune.cells.raw, subset = broad_label == "T cell")
primary.myeloid.cells.raw<-subset(primary.immune.cells.raw, subset = broad_label == "Myeloid")
primary.plasma.cells.raw<-subset(primary.immune.cells.raw, subset = broad_label == "Plasma cell")


normalize_raw<-function(data_subset){
  print(as.character(unique(data_subset$broad_label)))
  data_subset <- NormalizeData(data_subset)
  data_subset <- FindVariableFeatures(data_subset, selection.method = "vst", nfeatures = 2000)
  data_subset <- ScaleData(data_subset) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))
  
  # dimension reduction
  data_subset <- RunPCA(data_subset, ndims.print = 1:10, nfeatures.print = 10)
  data_subset <- RunUMAP(data_subset, dims = 1:30)
  data_subset <- RunTSNE(data_subset, dims = 1:30)
  
  ## cell cycle gene expression
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  data_subset <- CellCycleScoring(data_subset, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  
  ## regress out cell cycle and other covariates
  #Transformed data will be available in the SCT assay, which is set as the default after running sctransform
  #By default, sctransform accounts for cellular sequencing depth, or nUMIs.
  data_subset <- SCTransform(data_subset, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE)
  
  # dimension reduction
  data_subset <- RunPCA(data_subset, verbose = FALSE)
  data_subset <- RunUMAP(data_subset, dims = 1:30)
  data_subset <- RunTSNE(data_subset, dims = 1:30)
  
  # cluster
  data_subset <- FindNeighbors(data_subset, reduction = "pca", dims = 1:20)
  data_subset <- FindClusters(data_subset, resolution = 0.8)
  
  # relabel meta data
  data_subset@meta.data$orig.ident[which(data_subset@meta.data$orig.ident=="Ctrl")]<-"Control"
  data_subset@meta.data$orig.ident[which(data_subset@meta.data$orig.ident=="neo")]<-"Neonatal"
  
  print(as.character(unique(data_subset$broad_label)))
  print(table(data_subset@meta.data$elmentaite_annotation_V2,data_subset@meta.data$seurat_clusters))
  
  
  ## plot know cells on top of unknown (UMAP)
  umap_mat<-as.data.frame(Embeddings(object = data_subset, reduction = "umap"))#
  umap_mat$cell<-rownames(umap_mat)
  
  meta_primary<-data_subset@meta.data
  meta_primary$cell<-rownames(meta_primary)
  
  plt<-merge(umap_mat, meta_primary, by="cell")
  
  primary.immune.cells_known<-plt[which(plt$elmentaite_annotation_V2 != "unknown_filtered"),]
  primary.immune.cells_unknown<-plt[which(plt$elmentaite_annotation_V2 == "unknown_filtered"),]
  
  cell_labelled<-ggplot()+
    geom_point(aes(UMAP_1,UMAP_2, color=elmentaite_annotation_V2),primary.immune.cells_unknown,size=1.5)+
    geom_point(aes(UMAP_1,UMAP_2, color=elmentaite_annotation_V2),primary.immune.cells_known,size=1.5)+
    colcale_cols_manual+theme_classic()+
    th+theme(plot.title = element_text(size = 16, face = "bold"))
  
  clusters_labelled<-DimPlot(data_subset, reduction = "umap", pt.size=1, label=T)
  
  ggsave(file=paste(here("figs"),"/",as.character(unique(data_subset$broad_label)),"_UMAP_cluster.pdf",sep=""),
         grid.arrange(cell_labelled,clusters_labelled, ncol=2, widths=c(3,2)), w=15,h=5)
  ggsave(file=paste(here("figs/jpeg"),"/",as.character(unique(data_subset$broad_label)),"_UMAP_cluster.jpeg",sep=""),
         grid.arrange(cell_labelled,clusters_labelled, ncol=2, widths=c(3,2)), w=15,h=5)
  
  ggsave(file=paste(here("figs/jpeg"),"/",as.character(unique(data_subset$broad_label)),"_UMAP_cluster_quality.jpeg",sep=""),
  grid.arrange(FeaturePlot(data_subset, features = "nFeature_RNA",reduction = "umap", min.cutoff = "q9", pt.size=1),
               FeaturePlot(data_subset, features = "percent.mt",reduction = "umap", min.cutoff = "q9", pt.size=1),
               FeaturePlot(data_subset, features = "nCount_RNA",reduction = "umap", min.cutoff = "q9", pt.size=1)), w=5,h=15)
  
  DimPlot(data_subset, reduction = "umap", split.by = "orig.ident",label=T, pt.size=0.25)
  ggsave(file=paste(here("figs/jpeg"),"/",as.character(unique(data_subset$broad_label)),"_immune_primary_UMAP_source.jpeg",sep=""), w=10,h=4)
  
  FeaturePlot(data_subset, features = c("IGHA1","IGHG1","C1QA","ITGAX","SELL","TRDC","FCER2","POU2F2","CD3D"),reduction = "umap", min.cutoff = "q9", pt.size=1)
  ggsave(file=paste(here("figs/jpeg"),"/",as.character(unique(data_subset$broad_label)),"_immune_primary_UMAP_markers.jpeg",sep=""), w=22,h=20)
  
  FeaturePlot(data_subset, features = c("IL7R","CCR7","CD14","LYZ","IL7R","S100A4","MS4A1","CD8A","FCGR3A","MS4A7","GNLY","NKG7","FCER1A","CST3"),reduction = "umap", min.cutoff = "q9", pt.size=1)
  ggsave(file=paste(here("figs/jpeg"),"/",as.character(unique(data_subset$broad_label)),"_immune_primary_UMAP_seuratPBMC_markers.jpeg",sep=""), w=22,h=20)
  
  data_subset$pos_neg<-"whole"
  data_subset$pos_neg[grep("POS",data_subset$individual)]<-"POS"
  data_subset$pos_neg[grep("NEG",data_subset$individual)]<-"NEG"
  
  DimPlot(data_subset, reduction = "umap", split.by = "pos_neg",label=T, pt.size=0.25)
  ggsave(file=paste(here("figs/jpeg"),"/",as.character(unique(data_subset$broad_label)),"_immune_primary_UMAP_posneg.jpeg",sep=""), w=10,h=4)
  
  data_subset
  }


primary.B.cells<-normalize_raw(primary.B.cells.raw)
primary.T.cells<-normalize_raw(primary.T.cells.raw)
primary.myeloid.cells<-normalize_raw(primary.myeloid.cells.raw)
primary.plasma.cells<-normalize_raw(primary.plasma.cells.raw)



#save(primary.B.cells,primary.T.cells,primary.myeloid.cells,primary.plasma.cells, file = "../../data/d10x_immune_primary_normalized_seperately.RData")

primary.B.cells$cluster_ID<-"Unidentified"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==0)]<-"FCER2 B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==1)]<-"Memory B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==2)]<-"FCER2 B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==3)]<-"B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==4)]<-"FCER2 B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==5)]<-"FCER2 B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==6)]<-"Memory B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==7)]<-"FCER2 B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==8)]<-"Memory B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==9)]<-"Cycling B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==10)]<-"Memory B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==11)]<-"Cycling B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==12)]<-"Cycling B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==13)]<-"Cycling B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==14)]<-"B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==15)]<-"B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==16)]<-"CD4 T cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==17)]<-"Cycling B cell"
primary.B.cells$cluster_ID[which(primary.B.cells$seurat_clusters==18)]<-"Memory B cell"


primary.T.cells$cluster_ID<-"Unidentified"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==0)]<-"CD8 T cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==1)]<-"Tfh"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==2)]<-"CD4 T cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==3)]<-"CD4 T cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==4)]<-"Treg"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==5)]<-"Activated T"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==6)]<-"gd T/NK cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==7)]<-"CD4 T cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==8)]<-"Activated T"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==9)]<-"CD4 T cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==10)]<-"CD8 T cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==11)]<-"CD8 T cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==12)]<-"CD8 T cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==13)]<-"gd T/NK cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==14)]<-"gd T/NK cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==15)]<-"Cycling B cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==16)]<-"CD8 T cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==17)]<-"gd T/NK cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==18)]<-"gd T/NK cell"
primary.T.cells$cluster_ID[which(primary.T.cells$seurat_clusters==19)]<-"CD4 T cell"


primary.myeloid.cells$cluster_ID<-"Unidentified"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==0)]<-"Monocyte"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==1)]<-"cDC2"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==2)]<-"Macrophage"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==3)]<-"Macrophage"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==4)]<-"Macrophage"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==5)]<-"Monocyte"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==6)]<-"Monocyte"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==7)]<-"Monocyte"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==8)]<-"cDC2"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==9)]<-"cDC2"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==10)]<-"cDC1"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==11)]<-"CD4 T cell"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==12)]<-"cDC2"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==13)]<-"activated DC"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==14)]<-"pDC"
primary.myeloid.cells$cluster_ID[which(primary.myeloid.cells$seurat_clusters==15)]<-"Macrophage"



primary.plasma.cells$cluster_ID<-"Unidentified"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==0)]<-"IgA plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==1)]<-"IgA plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==2)]<-"Activated B cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==3)]<-"IgG plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==4)]<-"IgG plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==5)]<-"IgA plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==6)]<-"IgA plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==7)]<-"IgA plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==8)]<-"IgA plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==9)]<-"CD4 T cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==10)]<-"Activated B cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==11)]<-"IgG plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==12)]<-"IgG plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==13)]<-"IgA plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==14)]<-"IgA plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==15)]<-"IgG plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==16)]<-"IgG plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==17)]<-"IgG plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==18)]<-"IgG plasma cell"
primary.plasma.cells$cluster_ID[which(primary.plasma.cells$seurat_clusters==19)]<-"Venous endothelial cell"


immune_cell_labels<-rbind(primary.B.cells@meta.data[,c("seurat_clusters","cluster_ID")],
                          primary.T.cells@meta.data[,c("seurat_clusters","cluster_ID")],
                          primary.myeloid.cells@meta.data[,c("seurat_clusters","cluster_ID")],
                          primary.plasma.cells@meta.data[,c("seurat_clusters","cluster_ID")],
                          primary.mast.cells@meta.data[,c("seurat_clusters","cluster_ID")],
                          primary.S4.cells@meta.data[,c("seurat_clusters","cluster_ID")],
                          primary.myofibroblast.cells@meta.data[,c("seurat_clusters","cluster_ID")])




save(immune_cell_labels, file=here("output","immune_iterative_label.Rdata"))





#################
### cell label to entire UMAP
##################

## add cell label information
immune_cell_labels$index<-rownames(immune_cell_labels)
immune_cell_labels<-immune_cell_labels[match(colnames(primary.immune.cells), immune_cell_labels$index),]
identical(colnames(primary.immune.cells), immune_cell_labels$index)

primary.immune.cells<- AddMetaData(primary.immune.cells, immune_cell_labels$seurat_clusters, col.name = "seurat_interative_clustering")
primary.immune.cells<- AddMetaData(primary.immune.cells, immune_cell_labels$cluster_ID, col.name = "cell_type")

table(primary.immune.cells@meta.data$cell_type)

primary.immune.cells@meta.data$cell_type<-as.factor(primary.immune.cells@meta.data$cell_type)
primary.immune.cells@meta.data$cell_type<-factor(primary.immune.cells@meta.data$cell_type, levels = c(
  "Activated B cell","activated DC" ,"Activated T","Arterial endothelial cell" ,"B cell","CD4 T cell",
  "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
  "IgA plasma cell","IgG plasma cell","Lymphatic endothelial cell","Macrophage" ,"mast cells","Memory B cell" ,
  "Monocyte","pDC" ,"pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"TA","Tfh" ,
  "Treg", "Venous endothelial cell","Glial cell","myofibroblast","Mast cell","unknown_filtered","Doublet","Neonatal CD4 T cell","Neonatal B cell"))

#' ### dev cell labels
cols_manual<-c(  "#238443","#810f7c" ,"#02818a","#8c510a" ,"#78c679","#67a9cf",
                 "#3690c0","#8073ac","#b2abd2","#238443","#dd3497","#1d91c0","#addd8e","#014636",
                 "#253494","#081d58","#bf812d","#7a0177" ,"#ce1256","#006d2c" ,
                 "#a50f15","#542788" ,"#e31a1c","#e08214","#ef6548","#fdd49e" ,"#f46d43","#4575b4" ,
                 "#016c59", "#bf812d","#c51b7d","#fb9a99","#fb6a4a","lightgrey","black","#c6dbef","#c7e9c0")


names(cols_manual) <- levels(primary.immune.cells@meta.data$cell_type)
fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)




### PLot each subset with final labels
DimPlot(primary.B.cells, reduction = "umap", group.by = "cluster_ID", pt.size=1, label=T)+colcale_cols_manual
ggsave(file=here("figs","Bcell_primary_immune_UMAP_cluster_iterative_label.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","Bcell_primary_immune_UMAP_cluster_iterative_label.jpeg"), w=7,h=5)

DimPlot(primary.T.cells, reduction = "umap", group.by = "cluster_ID", pt.size=1, label=T)+colcale_cols_manual
ggsave(file=here("figs","Tcell_primary_immune_UMAP_cluster_iterative_label.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","Tcell_primary_immune_UMAP_cluster_iterative_label.jpeg"), w=7,h=5)

DimPlot(primary.myeloid.cells, reduction = "umap",  group.by = "cluster_ID",pt.size=1, label=T)+colcale_cols_manual
ggsave(file=here("figs","Myeloid_primary_immune_UMAP_cluster_iterative_label.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","Myeloid_primary_immune_UMAP_cluster_iterative_label.jpeg"), w=7,h=5)

DimPlot(primary.plasma.cells, reduction = "umap",  group.by = "cluster_ID",pt.size=1, label=T)+colcale_cols_manual
ggsave(file=here("figs","Plasma_primary_immune_UMAP_cluster_iterative_label.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","Plasma_primary_immune_UMAP_cluster_iterative_label.jpeg"), w=7,h=5)




## Plot all together
DimPlot(primary.immune.cells, reduction = "umap", group.by = "cell_type", pt.size=1, label=T)+colcale_cols_manual
ggsave(file=here("figs","primary_immune_UMAP_iterative_cell_type.pdf"), w=10,h=5)
ggsave(file=here("figs/jpeg","primary_immune_iterative_cell_type.jpeg"), w=10,h=5)

DimPlot(primary.immune.cells, reduction = "umap", group.by = "cell_type", pt.size=0.5)+colcale_cols_manual
ggsave(file=here("figs","primary_immune_UMAP_iterative_cell_type_no_label.pdf"), w=10,h=5)
ggsave(file=here("figs/jpeg","primary_immune_iterative_cell_type_no_label.jpeg"), w=10,h=5)