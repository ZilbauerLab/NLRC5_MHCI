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
identical(rownames(d10x.primary.raw@meta.data), cell_label$index)

d10x.primary.raw<- AddMetaData(d10x.primary.raw, cell_label$seurat_clusters, col.name = "seurat_wholePrimary_clusters")
d10x.primary.raw<- AddMetaData(d10x.primary.raw, cell_label$elmentaite_annotation_V2, col.name = "elmentaite_annotation_V2")


## subset to just epithelial
primary.epi.cells<-subset(d10x.primary.raw, subset = seurat_wholePrimary_clusters %in% c(12,14,25))

# An object of class Seurat 
# 24945 features across 2999 samples within 1 assay 
# Active assay: RNA (24945 features, 0 variable features)

## remove neonatal
primary.epi.cells<-subset(primary.epi.cells, subset = orig.ident != "neo")


primary.epi.cells <- NormalizeData(primary.epi.cells)
primary.epi.cells <- FindVariableFeatures(primary.epi.cells, selection.method = "vst", nfeatures = 2000)
primary.epi.cells <- ScaleData(primary.epi.cells) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
primary.epi.cells <- RunPCA(primary.epi.cells, ndims.print = 1:10, nfeatures.print = 10)
primary.epi.cells <- RunUMAP(primary.epi.cells, dims = 1:30)
primary.epi.cells <- RunTSNE(primary.epi.cells, dims = 1:30)



######################
## cell cycle gene expression
######################

## which PC associated to nfeature? 2
pca_mat<-as.data.frame(Embeddings(object = primary.epi.cells, reduction = "pca"))
pca_mat$cell<-rownames(pca_mat)
meta_primary<-primary.epi.cells@meta.data
meta_primary$cell<-rownames(meta_primary)
plt<-merge(pca_mat, meta_primary,by="cell")

sapply(2:51, function(x) cor(plt[,x],plt[,"nFeature_RNA"]))
FeaturePlot(primary.epi.cells, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
## PC2 tops genes include cell cycle genes: "TPX2"   "TOP2A"  "MKI67"  "CENPF"  "SMC4"   "TUBB4B"

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

primary.epi.cells <- CellCycleScoring(primary.epi.cells, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(primary.epi.cells, reduction="pca")
ggsave(file=here("figs","epi_primary_cell_cycle_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","epi_primary_cell_cycle_PCA.jpeg"), w=5,h=5)

FeaturePlot(primary.epi.cells, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","epi_primary_nfeature_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","epi_primary_nfeature_PCA.jpeg"), w=5.25,h=5)

DimPlot(primary.epi.cells, reduction="umap", group.by = "Phase")


## regress out cell cycle and other covariates
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#By default, sctransform accounts for cellular sequencing depth, or nUMIs.
primary.epi.cells <- SCTransform(primary.epi.cells, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE)


# dimension reduction
primary.epi.cells <- RunPCA(primary.epi.cells, verbose = FALSE)
primary.epi.cells <- RunUMAP(primary.epi.cells, dims = 1:30)
primary.epi.cells <- RunTSNE(primary.epi.cells, dims = 1:30)

# cluster
primary.epi.cells <- FindNeighbors(primary.epi.cells, reduction = "pca", dims = 1:20)
primary.epi.cells <- FindClusters(primary.epi.cells, resolution = 0.6)

# relabel meta data
primary.epi.cells@meta.data$orig.ident[which(primary.epi.cells@meta.data$orig.ident=="Ctrl")]<-"Control"
                                  #primary.epi.cells@meta.data$orig.ident[which(primary.epi.cells@meta.data$orig.ident=="neo")]<-"Neonatal"

saveRDS(primary.epi.cells, file = here("data","d10x_epi_primary_normalized.rds"))
# 
# 
# 
# 






primary.epi.cells<-readRDS(here("data","d10x_epi_primary_normalized.rds"))


## visualize
DimPlot(primary.epi.cells, reduction = "umap", group.by = "Phase", pt.size=1)
DimPlot(primary.epi.cells, reduction = "umap", group.by = "Phase",split.by = "orig.ident", pt.size=1)


DimPlot(primary.epi.cells, reduction = "umap", split.by = "orig.ident",label=T, pt.size=0.25)
ggsave(file=here("figs","epi_primary_UMAP_source.pdf"), w=20,h=4)
ggsave(file=here("figs/jpeg","epi_primary_UMAP_source.jpeg"), w=10,h=4)

DimPlot(primary.epi.cells, reduction = "tsne", split.by = "orig.ident",label=T, pt.size=0.25)
ggsave(file=here("figs","epi_primary_TSNE_raw_source.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","epi_primary_TSNE__source.jpeg"), w=7,h=6)

DimPlot(primary.epi.cells, reduction = "tsne",  pt.size=0.25)
ggsave(file=here("figs","epi_primary_TSNE_raw_cluster.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","epi_primary_TSNE_raw_cluster.jpeg"), w=7,h=6)

DimPlot(primary.epi.cells, reduction = "umap",  pt.size=0.25)
ggsave(file=here("figs","epi_primary_UMAP_raw_cluster.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","epi_primary_UMAPraw_cluster.jpeg"), w=7,h=6)


DimPlot(primary.epi.cells, reduction="pca", pt.size=1)

DimPlot(primary.epi.cells, reduction = "umap", group.by = "individual", pt.size=1)
DimPlot(primary.epi.cells, reduction = "umap", split.by = "individual", group.by = "orig.ident", pt.size=1)




table(primary.epi.cells@meta.data$elmentaite_annotation_V2)

primary.epi.cells@meta.data$elmentaite_annotation_V2<-as.factor(primary.epi.cells@meta.data$elmentaite_annotation_V2)
primary.epi.cells@meta.data$elmentaite_annotation_V2<-factor(primary.epi.cells@meta.data$elmentaite_annotation_V2, levels = c(
  "Arterial endothelial cell", "Venous endothelial cell",   
  "B cell","Memory B cell","Cycling myeloid cells","CD4 T cell","CD8 T cell","IgA plasma cell",       
  "Glial cell",   
  "BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine",
  "Goblet cell","IL2RG+ enterocyte (M cell)","Paneth cell","TA","Tuft","unknown_filtered"))
  
#' ### dev cell labels
cols_manual<-c( "#D3ADF3","#D3ADF3",
                "pink","pink","pink","pink","pink","pink",
                "#ADF3D2",
                "#3D0AF2","#6521C1","#67D364","#367C34",
                "#B921C1","#308AC8","#46B3DC",
                "#C86E30","#BFC121","#C12134","lightgrey")



DimPlot(primary.epi.cells, reduction = "umap", group.by = "elmentaite_annotation_V2", pt.size=1, label=T)+scale_color_manual(values=cols_manual)
ggsave(file=here("figs","primary_epithelial_UMAP_cellLabel_devcell.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","primary_epithelial_UMAP_cellLabel_devcell.jpeg"), w=7,h=5)

DimPlot(primary.epi.cells, reduction = "umap", group.by = "elmentaite_annotation_V2", pt.size=1, label=F)+scale_color_manual(values=cols_manual)
ggsave(file=here("figs","primary_epithelial_UMAP_cellLabel_devcell_nolab.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","primary_epithelial_UMAP_cellLabel_devcell_nolab.jpeg"), w=7,h=5)


DimPlot(primary.epi.cells, reduction = "umap", pt.size=1, label=T)
ggsave(file=here("figs","primary_epithelial_UMAP_cluster.pdf"), w=6,h=5)
ggsave(file=here("figs/jpeg","primary_epithelial_UMAP_cluster.jpeg"), w=6,h=5)


## relable clusters with cell types
table(primary.epi.cells@meta.data$elmentaite_annotation_V2,primary.epi.cells@meta.data$seurat_clusters)

## what type of t cells are mixed in 
FeaturePlot(primary.epi.cells, features = c("CD8A","IL7R","CCR7","S100A4","SELL"),reduction = "umap", min.cutoff = "q9", pt.size=1)
# expressing CD8A so CD8 T cell

### big dotplot
genes_df<-read.csv(here("data","single_cell_markers.csv"))
DotPlot(primary.epi.cells, features = unique(genes_df$Marker), group.by='seurat_clusters')+
  RotatedAxis()+theme(axis.title.x=element_blank())+ylab("seurat_clusters")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = unit(c(1,0.5,0.1,0.1), "cm"))+th_present

ggsave(file=here("figs","Allsamples_markerdotplot_celllabelled_primary_epithelial.pdf"), w=14,h=4)
ggsave(file=here("figs","Allsamples_markerdotplot_celllabelled_primary_epithelial.jpeg"), w=14,h=4)



primary.epi.cells$cluster_ID<-"Unidentified"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==0)]<-"Goblet cell"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==1)]<-"crypt"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==2)]<-"enterocyte"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==3)]<-"TA"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==4)]<-"Paneth"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==5)]<-"Goblet cell"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==6)]<-"CD8 T cell"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==7)]<-"Goblet cell"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==8)]<-"Paneth (UC only)"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==9)]<-"crypt"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==10)]<-"Memory B cell"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==11)]<-"enteroendocrine"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==12)]<-"enterocyte"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==13)]<-"Tuft"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==14)]<-"BEST4 enterocyte"
primary.epi.cells$cluster_ID[which(primary.epi.cells$seurat_clusters==15)]<-"enterocyte"





primary.epi.cells@meta.data$cluster_ID<-as.factor(primary.epi.cells@meta.data$cluster_ID)
primary.epi.cells@meta.data$cluster_ID<-factor(primary.epi.cells@meta.data$cluster_ID, levels = c(
  "Memory B cell", "CD8 T cell","Neonatal_cell","Unidentified_only_UC",
  "BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine",
  "Goblet cell", "Paneth","Paneth (UC only)","TA","Tuft"))





#' ### dev cell labels
cols_manual_less<-c( "pink","pink",
                "pink","pink",
                "#3D0AF2","#6521C1","#67D364","#367C34",
                "#B921C1","#308AC8",
                "#C86E30","#BFC121","#C12134")


DimPlot(primary.epi.cells, reduction = "umap", group.by = "cluster_ID", pt.size=1, label=T)
ggsave(file=here("figs","primary_epithelial_UMAP_label_expanded.pdf"), w=7,h=5)
ggsave(file=here("figs/jpeg","primary_epithelial_UMAP_label_expanded.jpeg"), w=7,h=5)


DimPlot(primary.epi.cells, reduction = "umap", split.by = "orig.ident",  group.by = "cluster_ID",label=F, pt.size=0.25)
ggsave(file=here("figs","epi_primary_UMAP_source_cluster.pdf"), w=10,h=4)
ggsave(file=here("figs/jpeg","epi_primary_UMAP_source_cluster.jpeg"), w=10,h=4)


save(primary.epi.cells, file=here("data","primary.epi.cells.RData"))

epi_cell_labels<-primary.epi.cells@meta.data[,c("seurat_clusters","cluster_ID")]
epi_cell_labels$seurat_clusters<-as.character(epi_cell_labels$seurat_clusters)
epi_cell_labels$cluster_ID<-as.character(epi_cell_labels$cluster_ID)

save(epi_cell_labels, file=here("output","epi_celltype_label.Rdata"))




