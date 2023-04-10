
#'---
#'title: Stimulation specific clustering
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
library(here)

options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))




d10x.stim.raw<-readRDS(here("data","d10x_raw_merged.rds"))

## subset no treatment
NT.cells<-subset(d10x.stim.raw, subset = orig.ident == "NT")

# An object of class Seurat 
# An object of class Seurat 
# 21548 features across 8649 samples within 1 assay 
# Active assay: RNA (21548 features, 0 variable features)

NT.cells <- NormalizeData(NT.cells)
NT.cells <- FindVariableFeatures(NT.cells, selection.method = "vst", nfeatures = 2000)
NT.cells <- ScaleData(NT.cells) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
NT.cells <- RunPCA(NT.cells, ndims.print = 1:10, nfeatures.print = 10)
NT.cells <- RunUMAP(NT.cells, dims = 1:30)
NT.cells <- RunTSNE(NT.cells, dims = 1:30)


######################
## cell cycle gene expression
######################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

NT.cells <- CellCycleScoring(NT.cells, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

NT.cells <- RunPCA(NT.cells)

DimPlot(NT.cells, reduction="pca")
ggsave(file=here("figs","NT_notintegrated_Phase_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","NT_notintegrated_Phase_PCA.jpeg"), w=5.25,h=5)
FeaturePlot(NT.cells, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","NT_notintegrated_nFeature_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","NT_notintegrated_nFeature_PCA.jpeg"), w=5.25,h=5)
FeaturePlot(NT.cells, features = "percent.mt",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","NT_notintegrated_perMT_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","NT_notintegrated_perMT_PCA.jpeg"), w=5.25,h=5)

PC_mat<-as.data.frame(Embeddings(object = NT.cells, reduction = "pca"))
PC_mat$cell<-rownames(PC_mat)

cor(NT.cells$nFeature_RNA, PC_mat$PC_2)
cor(NT.cells$nFeature_RNA, PC_mat$PC_1)
cor(NT.cells$percent.mt, PC_mat$PC_1)
cor(NT.cells$percent.mt, PC_mat$PC_2)

## going to regress phase and nFeature



##########################
## integrate across donor CCA with SCtransform on each donor
#########################
NT.cells.list <- SplitObject(NT.cells, split.by = "individual")
for (i in 1:length(NT.cells.list)) {
  suppressWarnings( NT.cells.list[[i]] <- SCTransform(NT.cells.list[[i]], vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE) )
}

#' Perform integration
NT.cells.features <- SelectIntegrationFeatures(object.list = NT.cells.list, nfeatures = 3000)

options(future.globals.maxSize = 3000 * 1024^2)
NT.cells.list <- PrepSCTIntegration(object.list = NT.cells.list, anchor.features = NT.cells.features,verbose = FALSE)

SCtransform.anchors <- FindIntegrationAnchors(object.list = NT.cells.list, normalization.method = "SCT",
                                              anchor.features = NT.cells.features, verbose = FALSE)
NT.cells.list.integrated <- IntegrateData(anchorset = SCtransform.anchors, normalization.method = "SCT",
                                          verbose = FALSE)


# # dimension reduction
NT.cells.list.integrated <- RunPCA(NT.cells.list.integrated, verbose = FALSE)
NT.cells.list.integrated <- RunUMAP(NT.cells.list.integrated, dims = 1:12)
NT.cells.list.integrated <- RunTSNE(NT.cells.list.integrated, dims = 1:12)

NT.cells.list.integrated

## Determine the right number of clusters (will also base on marker genes)
ElbowPlot(NT.cells.list.integrated)

                        # # After PCA, significant principal components were identified using
                        # # the permutation test 51 , implemented using the permutationPA function from the
                        # # jackstraw R package. This test identified 13 and 15 significant principal compo-
                        # #   nents in the 10X and SMART-Seq2 datasets
                        # NT.cells.list.integrated <- JackStraw(NT.cells.list.integrated, num.replicate = 100)
                        # NT.cells.list.integrated <- ScoreJackStraw(NT.cells.list.integrated, dims = 1:20)
                        # JackStrawPlot(NT.cells.list.integrated, dims = 1:20)
                        # ## seems like a drop in significance after PC13


# cluster
NT.cells.list.integrated <- FindNeighbors(NT.cells.list.integrated, reduction = "pca", dims = 1:13)
NT.cells.list.integrated <- FindClusters(NT.cells.list.integrated, resolution = 0.25)


saveRDS(NT.cells.list.integrated, file = here("data","NT_normalized_integrated.rds"))





#######################################################################################################################################################

NT.cells.integrated<-readRDS(here("data","NT_normalized_integrated.rds"))
NT.cells.integrated





DimPlot(NT.cells.integrated, reduction = "tsne", pt.size=1, label=T)
DimPlot(NT.cells.integrated, reduction = "umap", pt.size=1, label=T)


## visualize
DimPlot(NT.cells.integrated, reduction = "umap", pt.size=1, label=T)
ggsave(file=here("figs","UMAP_NT_integrated_clusters.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_NT_integrated_clusters.jpeg"), w=5,h=4)


DimPlot(NT.cells.integrated, reduction = "umap", group.by = "Phase", pt.size=1)
ggsave(file=here("figs","UMAP_NT_integrated_phase.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_NT_integrated_phase.jpeg"), w=5,h=4)

DimPlot(NT.cells.integrated, reduction="pca", pt.size=1)
ggsave(file=here("figs","PCA_NT_integrated.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_NT_integrated.jpeg"), w=5,h=4)

DimPlot(NT.cells.integrated, reduction = "tsne", label = T, pt.size=0.5)
ggsave(file=here("figs","tsne_NT_integrated.pdf"), w=7,h=6)
ggsave(file=here("figs/jpeg","tsne_NT_integrated.jpeg"), w=7,h=6)

DimPlot(NT.cells.integrated, reduction = "umap", group.by = "individual", pt.size=1)
ggsave(file=here("figs","UMAP_NT_integrated_individual.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_NT_integrated_individual.jpeg"), w=5,h=4)

DimPlot(NT.cells.integrated, reduction = "umap", split.by = "individual",  pt.size=1)
ggsave(file=here("figs","UMAP_NT_integrated_individualsplit.pdf"), w=15,h=4)
ggsave(file=here("figs/jpeg","UMAP_NT_integrated_individualsplit.jpeg"), w=15,h=4)


# ## Low quality cluster?
FeaturePlot(NT.cells.integrated, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","UMAP_NT_integrated_nFeature.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_NT_integrated_nFeature.jpeg"), w=5,h=4)

FeaturePlot(NT.cells.integrated, features = "percent.mt", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","UMAP_NT_integrated_MT.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_NT_integrated_MT.jpeg"), w=5,h=4)


FeaturePlot(NT.cells.integrated, features = "nFeature_RNA", reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","PCA_NT_integrated_nFeature.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_NT_integrated_nFeature.jpeg"), w=5,h=4)

FeaturePlot(NT.cells.integrated, features = "percent.mt", reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","PCA_NT_integrated_percentmt.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_NT_integrated_percentmt.jpeg"), w=5,h=4)

DimPlot(NT.cells.integrated, reduction = "pca", group.by = "Phase", pt.size=1)
ggsave(file=here("figs","PCA_NT_integrated_phase.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_NT_integrated_phase.jpeg"), w=5,h=4)

## genes in first pcs
pc_load<-Loadings(NT.cells.integrated[['pca']])
print(NT.cells.integrated[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(NT.cells.integrated, dims = 1:2, reduction = "pca")





##########################################################
## IFNg
##########################################################
## subset IFNg
IFNg.cells<-subset(d10x.stim.raw, subset = orig.ident == "IFNg")

# An object of class Seurat
# An object of class Seurat
# 21548 features across 8649 samples within 1 assay
# Active assay: RNA (21548 features, 0 variable features)

IFNg.cells <- NormalizeData(IFNg.cells)
IFNg.cells <- FindVariableFeatures(IFNg.cells, selection.method = "vst", nfeatures = 2000)
IFNg.cells <- ScaleData(IFNg.cells) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
IFNg.cells <- RunPCA(IFNg.cells, ndims.print = 1:10, nfeatures.print = 10)
IFNg.cells <- RunUMAP(IFNg.cells, dims = 1:30)
IFNg.cells <- RunTSNE(IFNg.cells, dims = 1:30)


######################
## cell cycle gene expression
######################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

IFNg.cells <- CellCycleScoring(IFNg.cells, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

IFNg.cells <- RunPCA(IFNg.cells)

DimPlot(IFNg.cells, reduction="pca")
ggsave(file=here("figs","IFNg_notintegrated_Phase_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","IFNg_notintegrated_Phase_PCA.jpeg"), w=5.25,h=5)
FeaturePlot(IFNg.cells, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","IFNg_notintegrated_nFeature_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","IFNg_notintegrated_nFeature_PCA.jpeg"), w=5.25,h=5)
FeaturePlot(IFNg.cells, features = "percent.mt",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","IFNg_notintegrated_perMT_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","IFNg_notintegrated_perMT_PCA.jpeg"), w=5.25,h=5)

PC_mat<-as.data.frame(Embeddings(object = IFNg.cells, reduction = "pca"))
PC_mat$cell<-rownames(PC_mat)

cor(IFNg.cells$nFeature_RNA, PC_mat$PC_2)
cor(IFNg.cells$nFeature_RNA, PC_mat$PC_1)
cor(IFNg.cells$percent.mt, PC_mat$PC_1)
cor(IFNg.cells$percent.mt, PC_mat$PC_2)

## going to regress phase and nFeature


##########################
## integrate across donor CCA with SCtransform on each donor
#########################
IFNg.cells.list <- SplitObject(IFNg.cells, split.by = "individual")
for (i in 1:length(IFNg.cells.list)) {
  suppressWarnings( IFNg.cells.list[[i]] <- SCTransform(IFNg.cells.list[[i]], vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE) )
}

#' Perform integration
IFNg.cells.features <- SelectIntegrationFeatures(object.list = IFNg.cells.list, nfeatures = 3000)

options(future.globals.maxSize = 3000 * 1024^2)
IFNg.cells.list <- PrepSCTIntegration(object.list = IFNg.cells.list, anchor.features = IFNg.cells.features,verbose = FALSE)

SCtransform.anchors <- FindIntegrationAnchors(object.list = IFNg.cells.list, normalization.method = "SCT",
                                              anchor.features = IFNg.cells.features, verbose = FALSE)
IFNg.cells.list.integrated <- IntegrateData(anchorset = SCtransform.anchors, normalization.method = "SCT",
                                            verbose = FALSE)

# # dimension reduction
IFNg.cells.list.integrated <- RunPCA(IFNg.cells.list.integrated, verbose = FALSE)
IFNg.cells.list.integrated <- RunUMAP(IFNg.cells.list.integrated, dims = 1:12)
IFNg.cells.list.integrated <- RunTSNE(IFNg.cells.list.integrated, dims = 1:12)

## Cluster (But will want to project on NT clusters, can then compared to unsupervised clusters)
IFNg.cells.list.integrated <- FindNeighbors(IFNg.cells.list.integrated, reduction = "pca", dims = 1:12)
IFNg.cells.list.integrated <- FindClusters(IFNg.cells.list.integrated, resolution = 0.25)

saveRDS(IFNg.cells.list.integrated, file = here("data","IFNg_normalized_integrated.rds"))






IFNg.cells.integrated<-readRDS(here("data","IFNg_normalized_integrated.rds"))
IFNg.cells.integrated


DimPlot(IFNg.cells.integrated, reduction = "umap", pt.size=1, label=T)
DimPlot(IFNg.cells.integrated, reduction = "tsne", pt.size=1, label=T)


## visualize
DimPlot(IFNg.cells.integrated, reduction = "umap", pt.size=1, label=T)
ggsave(file=here("figs","UMAP_IFNg_integrated_clusters.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_IFNg_integrated_clusters.jpeg"), w=5,h=4)


DimPlot(IFNg.cells.integrated, reduction = "umap", group.by = "Phase", pt.size=1)
ggsave(file=here("figs","UMAP_IFNg_integrated_phase.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_IFNg_integrated_phase.jpeg"), w=5,h=4)

DimPlot(IFNg.cells.integrated, reduction="pca", pt.size=1)
ggsave(file=here("figs","PCA_IFNg_integrated.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_IFNg_integrated.jpeg"), w=5,h=4)

DimPlot(IFNg.cells.integrated, reduction = "umap", group.by = "individual", pt.size=1)
ggsave(file=here("figs","UMAP_IFNg_integrated_individual.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_IFNg_integrated_individual.jpeg"), w=5,h=4)

DimPlot(IFNg.cells.integrated, reduction = "umap", split.by = "individual",  pt.size=1)
ggsave(file=here("figs","UMAP_IFNg_integrated_individualsplit.pdf"), w=15,h=4)
ggsave(file=here("figs/jpeg","UMAP_IFNg_integrated_individualsplit.jpeg"), w=15,h=4)


# ## Low quality cluster?
FeaturePlot(IFNg.cells.integrated, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","UMAP_IFNg_integrated_nFeature.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_IFNg_integrated_nFeature.jpeg"), w=5,h=4)

FeaturePlot(IFNg.cells.integrated, features = "percent.mt", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","UMAP_IFNg_integrated_MT.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_IFNg_integrated_MT.jpeg"), w=5,h=4)
#
FeaturePlot(IFNg.cells.integrated, features = "nFeature_RNA", reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","PCA_IFNg_integrated_nFeature.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_IFNg_integrated_nFeature.jpeg"), w=5,h=4)

FeaturePlot(IFNg.cells.integrated, features = "percent.mt", reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","PCA_IFNg_integrated_percentmt.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_IFNg_integrated_percentmt.jpeg"), w=5,h=4)

DimPlot(IFNg.cells.integrated, reduction = "pca", group.by = "Phase", pt.size=1)
ggsave(file=here("figs","PCA_IFNg_integrated_phase.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_IFNg_integrated_phase.jpeg"), w=5,h=4)






##########################################################
## TNFa
##########################################################
## subset TNFa
TNFa.cells<-subset(d10x.stim.raw, subset = orig.ident == "TNFa")

# An object of class Seurat
# An object of class Seurat
# 21548 features across 8649 samples within 1 assay
# Active assay: RNA (21548 features, 0 variable features)

TNFa.cells <- NormalizeData(TNFa.cells)
TNFa.cells <- FindVariableFeatures(TNFa.cells, selection.method = "vst", nfeatures = 2000)
TNFa.cells <- ScaleData(TNFa.cells) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
TNFa.cells <- RunPCA(TNFa.cells, ndims.print = 1:10, nfeatures.print = 10)
TNFa.cells <- RunUMAP(TNFa.cells, dims = 1:30)
TNFa.cells <- RunTSNE(TNFa.cells, dims = 1:30)


######################
## cell cycle gene expression
######################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

TNFa.cells <- CellCycleScoring(TNFa.cells, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

TNFa.cells <- RunPCA(TNFa.cells)

DimPlot(TNFa.cells, reduction="pca")
ggsave(file=here("figs","TNFa_notintegrated_Phase_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","TNFa_notintegrated_Phase_PCA.jpeg"), w=5.25,h=5)
FeaturePlot(TNFa.cells, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","TNFa_notintegrated_nFeature_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","TNFa_notintegrated_nFeature_PCA.jpeg"), w=5.25,h=5)
FeaturePlot(TNFa.cells, features = "percent.mt",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","TNFa_notintegrated_perMT_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","TNFa_notintegrated_perMT_PCA.jpeg"), w=5.25,h=5)

PC_mat<-as.data.frame(Embeddings(object = TNFa.cells, reduction = "pca"))
PC_mat$cell<-rownames(PC_mat)

cor(TNFa.cells$nFeature_RNA, PC_mat$PC_2)
cor(TNFa.cells$nFeature_RNA, PC_mat$PC_1)
cor(TNFa.cells$percent.mt, PC_mat$PC_1)
cor(TNFa.cells$percent.mt, PC_mat$PC_2)

## going to regress phase and nFeature


##########################
## integrate across donor CCA with SCtransform on each donor
#########################
TNFa.cells.list <- SplitObject(TNFa.cells, split.by = "individual")
for (i in 1:length(TNFa.cells.list)) {
  suppressWarnings( TNFa.cells.list[[i]] <- SCTransform(TNFa.cells.list[[i]], vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE) )
}

#' Perform integration
TNFa.cells.features <- SelectIntegrationFeatures(object.list = TNFa.cells.list, nfeatures = 3000)

options(future.globals.maxSize = 3000 * 1024^2)
TNFa.cells.list <- PrepSCTIntegration(object.list = TNFa.cells.list, anchor.features = TNFa.cells.features,verbose = FALSE)

SCtransform.anchors <- FindIntegrationAnchors(object.list = TNFa.cells.list, normalization.method = "SCT",
                                              anchor.features = TNFa.cells.features, verbose = FALSE)
TNFa.cells.list.integrated <- IntegrateData(anchorset = SCtransform.anchors, normalization.method = "SCT",
                                            verbose = FALSE)


# # dimension reduction
TNFa.cells.list.integrated <- RunPCA(TNFa.cells.list.integrated, verbose = FALSE)
TNFa.cells.list.integrated <- RunUMAP(TNFa.cells.list.integrated, dims = 1:12)
TNFa.cells.list.integrated <- RunTSNE(TNFa.cells.list.integrated, dims = 1:12)

## Cluster (But will want to project on NT clusters, can then compared to unsupervised clusters)
TNFa.cells.list.integrated <- FindNeighbors(TNFa.cells.list.integrated, reduction = "pca", dims = 1:12)
TNFa.cells.list.integrated <- FindClusters(TNFa.cells.list.integrated, resolution = 0.25)

saveRDS(TNFa.cells.list.integrated, file = here("data","TNFa_normalized_integrated.rds"))


TNFa.cells.integrated<-readRDS(here("data","TNFa_normalized_integrated.rds"))
TNFa.cells.integrated


DimPlot(TNFa.cells.integrated, reduction = "umap", pt.size=1, label=T)
DimPlot(TNFa.cells.integrated, reduction = "tsne", pt.size=1, label=T)



## visualize
DimPlot(TNFa.cells.integrated, reduction = "umap", pt.size=1, label=T)
ggsave(file=here("figs","UMAP_TNFa_integrated_clusters.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_TNFa_integrated_clusters.jpeg"), w=5,h=4)


DimPlot(TNFa.cells.integrated, reduction = "umap", group.by = "Phase", pt.size=1)
ggsave(file=here("figs","UMAP_TNFa_integrated_phase.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_TNFa_integrated_phase.jpeg"), w=5,h=4)

DimPlot(TNFa.cells.integrated, reduction="pca", pt.size=1)
ggsave(file=here("figs","PCA_TNFa_integrated.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_TNFa_integrated.jpeg"), w=5,h=4)

DimPlot(TNFa.cells.integrated, reduction = "umap", group.by = "individual", pt.size=1)
ggsave(file=here("figs","UMAP_TNFa_integrated_individual.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_TNFa_integrated_individual.jpeg"), w=5,h=4)

DimPlot(TNFa.cells.integrated, reduction = "umap", split.by = "individual",  pt.size=1)
ggsave(file=here("figs","UMAP_TNFa_integrated_individualsplit.pdf"), w=15,h=4)
ggsave(file=here("figs/jpeg","UMAP_TNFa_integrated_individualsplit.jpeg"), w=15,h=4)


# ## Low quality cluster?
FeaturePlot(TNFa.cells.integrated, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","UMAP_TNFa_integrated_nFeature.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_TNFa_integrated_nFeature.jpeg"), w=5,h=4)

FeaturePlot(TNFa.cells.integrated, features = "percent.mt", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","UMAP_TNFa_integrated_MT.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_TNFa_integrated_MT.jpeg"), w=5,h=4)

FeaturePlot(TNFa.cells.integrated, features = c("LGR5","EPCAM","FABP1"), min.cutoff = "q9")


FeaturePlot(TNFa.cells.integrated, features = "nFeature_RNA", reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","PCA_TNFa_integrated_nFeature.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_TNFa_integrated_nFeature.jpeg"), w=5,h=4)

FeaturePlot(TNFa.cells.integrated, features = "percent.mt", reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","PCA_TNFa_integrated_percentmt.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_TNFa_integrated_percentmt.jpeg"), w=5,h=4)

DimPlot(TNFa.cells.integrated, reduction = "pca", group.by = "Phase", pt.size=1)
ggsave(file=here("figs","PCA_TNFa_integrated_phase.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","PCA_TNFa_integrated_phase.jpeg"), w=5,h=4)




#######################
## presentation dotplot
#######################
genes_df<-data.frame(gene=c("LGR5","ASCL2","FABP1","KRT19","HELLS","PCNA","MKI67","TOP2A", "FCGBP","SPINK4", "SOX4","CD24"),
                     cellType=c("Stem","Stem","Early E","Early E","TA1","TA1","TA2","TA2","Goblet","Goblet","Tuft","Tuft"))

rect<-ggplot()+geom_rect(aes(xmin=1:nrow(genes_df), xmax=1:nrow(genes_df)+1, ymin=1, ymax=2,
                             fill=cellType),genes_df, color="black") + geom_hline(yintercept=-1, color="white")+
  theme_bw()+theme(legend.position = "none",
                   axis.title=element_blank(),
                   axis.text=element_blank(),
                   axis.ticks=element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_blank(),
                   panel.border = element_blank(),plot.margin = unit(c(0,5,0,0.75), "cm"))+
  geom_text(aes(x=seq(1,nrow(genes_df), 2)+1, y=0, label=genes_df$cellType[seq(1,nrow(genes_df), 2)]), size=4)+
  scale_fill_manual(values=c("#cf92beff","#364557","#9dbc26ff","#efe532ff","#F6A317","#89DFE2"))#



ggsave(grid.arrange(
  DotPlot(NT.cells.integrated, features = genes_df$gene)+
    RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                            plot.margin = unit(c(0.3,0.5,0.1,0.1), "cm")),
  rect,
  heights=c(0.9,0.1)),
  file=here("figs/jpeg","NT_markerdotplot.jpeg"), w=7,h=4)

ggsave(grid.arrange(
  DotPlot(NT.cells.integrated, features = genes_df$gene)+
    RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                            plot.margin = unit(c(0.3,0.5,0.1,0.1), "cm")),
  rect,
  heights=c(0.9,0.1)),
  file=here("figs","NT_markerdotplot.pdf"), w=7,h=4)





#'## R Session Info
sessionInfo()
