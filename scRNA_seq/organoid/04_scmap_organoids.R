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
library(scmap)
library(SingleCellExperiment)
library(scater)

options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))


##############################
## control organoids
##############################
NT.cells.integrated<-readRDS(here("data","NT_normalized_integrated.rds"))
          # NT.cells.integrated <- FindNeighbors(NT.cells.integrated, reduction = "pca", dims = 1:13)
          # NT.cells.integrated <- FindClusters(NT.cells.integrated, resolution = 0.175)
DefaultAssay(object = NT.cells.integrated) <- "SCT"


NT.cells.sce <- as.SingleCellExperiment(NT.cells.integrated)
counts(NT.cells.sce)<-as.matrix(counts(NT.cells.sce))

p1 <- plotExpression(NT.cells.sce, features = "HLA-E", x = "seurat_clusters") + theme(axis.text.x = element_text(angle = 45, 
                                                                                                   hjust = 1))
p2 <- plotPCA(NT.cells.sce, colour_by = "ident")
p1 + p2


logcounts(NT.cells.sce) <- log2(counts(NT.cells.sce) + 1)
# use gene names as feature symbols
rowData(NT.cells.sce)$feature_symbol <- rownames(NT.cells.sce)
# remove features with duplicated names
NT.cells.sce <- NT.cells.sce[!duplicated(rownames(NT.cells.sce)), ]
NT.cells.sce
NT.cells.sce <- selectFeatures(NT.cells.sce, suppress_plot = FALSE)

NT.cells.sce <- indexCluster(NT.cells.sce, "seurat_clusters")
heatmap(as.matrix(metadata(NT.cells.sce)$scmap_cluster_index))




##############################
## primary epithelial
##############################
load(here("data","primary.epi.cells.RData"))
primary.epi.cells_subset<-subset(primary.epi.cells, subset = individual %in% c(
  "T017","T019","T176","T189","T197","T202","T203","T024","T036","T44","T057",
  "T160","T161","T175","T182","T184","T180"))
rm(primary.epi.cells)

DefaultAssay(object = primary.epi.cells_subset) <- "SCT"


primary.epi.cells_subset.sce <- as.SingleCellExperiment(primary.epi.cells_subset)
counts(primary.epi.cells_subset.sce)<-as.matrix(counts(primary.epi.cells_subset.sce))

p1 <- plotExpression(primary.epi.cells_subset.sce, features = "HLA-E", x = "cluster_ID") + theme(axis.text.x = element_text(angle = 45, 
                                                                                                                 hjust = 1))
p2 <- plotPCA(primary.epi.cells_subset.sce, colour_by = "ident")
p1 + p2


logcounts(primary.epi.cells_subset.sce) <- log2(counts(primary.epi.cells_subset.sce) + 1)
# use gene names as feature symbols
rowData(primary.epi.cells_subset.sce)$feature_symbol <- rownames(primary.epi.cells_subset.sce)
# remove features with duplicated names
primary.epi.cells_subset.sce <- primary.epi.cells_subset.sce[!duplicated(rownames(primary.epi.cells_subset.sce)), ]
primary.epi.cells_subset.sce
primary.epi.cells_subset.sce <- selectFeatures(primary.epi.cells_subset.sce, suppress_plot = FALSE)

# ## cell markers only
 # rowData(primary.epi.cells_subset.sce)$scmap_features <- (rowData(primary.epi.cells_subset.sce)$feature_symbol%in%c("LGR5","ASCL2","FABP1","KRT19","HELLS","PCNA","MKI67","TOP2A", "FCGBP","SPINK4", "SOX4","CD24"))
 # table(rowData(primary.epi.cells_subset.sce)$scmap_features)


primary.epi.cells_subset.sce <- indexCluster(primary.epi.cells_subset.sce, "cluster_ID")
heatmap(as.matrix(metadata(primary.epi.cells_subset.sce)$scmap_cluster_index))


##############################
## Fujii Organoids
##############################
counts_matrix_filename = paste0(here("data"),"GSM3389580_SB_LogNormalized.txt")
counts_matrix <- read.table(file =  counts_matrix_filename, header = TRUE, row.names = 1, sep = "\t", as.is = TRUE)

counts <- CreateSeuratObject(counts = counts_matrix, project = "fujii")  # Seurat function to read in 10x count data
counts <- FindVariableFeatures(counts, selection.method = "vst")
counts <- ScaleData(counts) 
counts <- RunPCA(counts, verbose = FALSE)
counts <- RunUMAP(counts, dims = 1:25)
counts <- RunTSNE(counts, dims = 1:25)
counts <- FindNeighbors(counts, reduction = "pca", dims = 1:25)
counts <- FindClusters(counts, resolution = 0.5)

counts$cluster_ID<-"Unidentified"
counts$cluster_ID[which(counts$seurat_clusters==0)]<-"Stem"
counts$cluster_ID[which(counts$seurat_clusters==1)]<-"Early E"
counts$cluster_ID[which(counts$seurat_clusters==2)]<-"TA1"
counts$cluster_ID[which(counts$seurat_clusters==3)]<-"TA1"
counts$cluster_ID[which(counts$seurat_clusters==4)]<-"TA2"
counts$cluster_ID[which(counts$seurat_clusters==5)]<-"TA2"
counts$cluster_ID[which(counts$seurat_clusters==6)]<-"Goblet"


fujii.sce <- as.SingleCellExperiment(counts)
counts(fujii.sce)<-as.matrix(counts(fujii.sce))

p1 <- plotExpression(fujii.sce, features = "HLA-E", x = "cluster_ID") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p2 <- plotPCA(fujii.sce, colour_by = "ident")
p1 + p2


logcounts(fujii.sce) <- log2(counts(fujii.sce) + 1)
# use gene names as feature symbols
rowData(fujii.sce)$feature_symbol <- rownames(fujii.sce)
# remove features with duplicated names
fujii.sce <- fujii.sce[!duplicated(rownames(fujii.sce)), ]
fujii.sce
fujii.sce <- selectFeatures(fujii.sce, suppress_plot = FALSE)
table(rowData(fujii.sce)$scmap_features)

# ## cell markers only
# rowData(fujii.sce)$scmap_features <- (rowData(fujii.sce)$feature_symbol%in%c("LGR5","ASCL2","FABP1","KRT19","HELLS","PCNA","MKI67","TOP2A", "FCGBP","SPINK4", "SOX4","CD24"))
# table(rowData(fujii.sce)$scmap_features)


fujii.sce <- indexCluster(fujii.sce, "cluster_ID")
heatmap(as.matrix(metadata(fujii.sce)$scmap_cluster_index))



##############################
#'## scmap at cluster level
##############################
# primary_to_NT <- scmapCluster(
#   projection = primary.epi.cells_subset.sce,
#   index_list = list(
#     NT.cells.sce = metadata(NT.cells.sce)$scmap_cluster_index
#   )
# )

primary_to_NT <- scmapCluster(
  projection = NT.cells.sce,
  index_list = list(
    primary.epi.cells_subset.sce = metadata(primary.epi.cells_subset.sce)$scmap_cluster_index
  )
)

table(colData(NT.cells.sce)$seurat_clusters, primary_to_NT$scmap_cluster_labs)
plot(getSankey(colData(NT.cells.sce)$seurat_clusters,  primary_to_NT$scmap_cluster_labs[,1], plot_height=400))



fujii_to_NT <- scmapCluster(
  projection = NT.cells.sce,
  index_list = list(
    fujii.sce = metadata(fujii.sce)$scmap_cluster_index
  )
)

table(colData(NT.cells.sce)$seurat_clusters, fujii_to_NT$scmap_cluster_labs)
plot(getSankey(colData(NT.cells.sce)$seurat_clusters,  fujii_to_NT$scmap_cluster_labs[,1], plot_height=400))

head(fujii_to_NT$scmap_cluster_labs)
head(fujii_to_NT$scmap_cluster_siml)

##################################
## scmap individual cell
##################################
set.seed(1)
fujii.sce <- indexCell(fujii.sce)
names(metadata(fujii.sce)$scmap_cell_index)
length(metadata(fujii.sce)$scmap_cell_index$subcentroids)
dim(metadata(fujii.sce)$scmap_cell_index$subcentroids[[1]])
metadata(fujii.sce)$scmap_cell_index$subcentroids[[1]][,1:5]
dim(metadata(fujii.sce)$scmap_cell_index$subclusters)
metadata(fujii.sce)$scmap_cell_index$subclusters[1:5,1:5]

scmapCell_results <- scmapCell(
  NT.cells.sce, 
  list(
    fujii.sce = metadata(fujii.sce)$scmap_cell_index
  ))

names(scmapCell_results)

scmapCell_results$fujii.sce$cells[,1:3]

scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, 
  list(
    as.character(colData(fujii.sce)$cluster_ID)
  ))

head(scmapCell_clusters$scmap_cluster_labs)
plot(
  getSankey(
    colData(NT.cells.sce)$seurat_clusters, 
    scmapCell_clusters$scmap_cluster_labs[,"fujii.sce"],
    plot_height = 400
  )
)
