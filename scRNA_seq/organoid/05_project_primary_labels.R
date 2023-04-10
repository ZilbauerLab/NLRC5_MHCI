
#'---
#'title: Labelling Organoid Cell Clusters
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
#library(DirichletReg)



options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))


## primary epithelial
load(here("data","primary.epi.cells.RData"))
primary.epi.cells_subset<-subset(primary.epi.cells, subset = individual %in% c(
  "T017","T019","T176","T189","T197","T202","T203","T024","T036","T44","T057",
  "T160","T161","T175","T182","T184","T180"))
rm(primary.epi.cells)

primary.epi.cells_subset<-subset(primary.epi.cells_subset, subset = cluster_ID %in% c("BEST4 enterocyte","crypt","enterocyte","enteroendocrine",
                                                                                   "Goblet cell","Paneth","Paneth (UC only)",
                                                                                   "TA","Tuft" ))

##################################
#' # control organoids
##################################
NT.cells.integrated<-readRDS(here("data","NT_normalized_integrated.rds"))
# NT.cells.integrated <- FindNeighbors(NT.cells.integrated, reduction = "pca", dims = 1:13)
# NT.cells.integrated <- FindClusters(NT.cells.integrated, resolution = 0.175)

print(table(NT.cells.integrated@meta.data$seurat_clusters))
primary.epi.cells_subset@meta.data$cluster_ID<-as.character(primary.epi.cells_subset@meta.data$cluster_ID)
print(table(primary.epi.cells_subset@meta.data$cluster_ID))


#####################
#' ## Seurat label transfer
#####################
org_pri.anchors <- FindTransferAnchors(reference = primary.epi.cells_subset, query = NT.cells.integrated, dims = 1:20, project.query = TRUE)
predictions_primary <- TransferData(anchorset = org_pri.anchors, refdata = primary.epi.cells_subset$cluster_ID, dims = 1:20, k.weight=10)
predictions_primary
org_pri.cells.integrated <- AddMetaData(NT.cells.integrated, metadata = predictions_primary)



grid.arrange(DimPlot(org_pri.cells.integrated, reduction = "umap", pt.size=1,  label=F),DimPlot(org_pri.cells.integrated, reduction = "umap", pt.size=1, group.by = "predicted.id",  label=F), ncol=2)
ggsave(file=here("figs","UMAP_NT_labeltransfer_clusters.pdf"),grid.arrange(DimPlot(org_pri.cells.integrated, reduction = "umap", pt.size=0.5,  label=F),DimPlot(org_pri.cells.integrated, reduction = "umap", pt.size=0.5, group.by = "predicted.id",  label=F), ncol=2, widths=c(0.8,1)), w=11,h=4)
ggsave(file=here("figs/jpeg","UMAP_NT_labeltransfer_clusters.jpeg"),grid.arrange(DimPlot(org_pri.cells.integrated, reduction = "umap", pt.size=0.25,  label=F),DimPlot(org_pri.cells.integrated, reduction = "umap", pt.size=0.25, group.by = "predicted.id",  label=F), ncol=2, widths=c(0.8,1)), w=11,h=4)


grid.arrange(DimPlot(org_pri.cells.integrated, reduction = "tsne", pt.size=1,  label=F),
             DimPlot(org_pri.cells.integrated, reduction = "tsne", pt.size=1, group.by = "predicted.id",  label=F), ncol=2, widths=c(1,1.3))


table(org_pri.cells.integrated$seurat_clusters, org_pri.cells.integrated$predicted.id)

org_pri.cells.integrated@meta.data$predictions_primary<-org_pri.cells.integrated@meta.data$predicted.id

pri_pred<-org_pri.cells.integrated@meta.data[,c("seurat_clusters","predictions_primary")]
save(pri_pred, file=here("output","primary_predict_transfer.RData"))


#'  Bar plot
cnt_plt<-as.data.frame(table(org_pri.cells.integrated$seurat_clusters, org_pri.cells.integrated$predicted.id))

ggplot(cnt_plt, aes(x=Var1, y=Freq,fill=Var2 )) + geom_bar(position="stack", stat="identity", color="black")+theme_bw()+th+xlab("Cluster")+ylab("Cell Number")+
  scale_fill_manual(values=c("#addd8e","#e31a1c","#2171b5","#238443","#74c476","#addd8e","#6baed6","#ffff33","#00441b"), name="Primary Epi Type")
ggsave(file=here("figs","NT_seurat_transfer_clusters_bar.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","NT_seurat_transfer_clusters_bar.jpeg"),w=6,h=4)



################
#'## Integrating organoid and primary and clustering together
################
combined_primary_NT<-list(NT.cells.integrated, primary.epi.cells_subset)

combined.features <- SelectIntegrationFeatures(object.list = combined_primary_NT, nfeatures = 3000)

options(future.globals.maxSize = 3000 * 1024^2)
combined_primary_NT <- PrepSCTIntegration(object.list = combined_primary_NT, anchor.features = combined.features,verbose = FALSE)

SCtransform.anchors <- FindIntegrationAnchors(object.list = combined_primary_NT, normalization.method = "SCT",
                                              anchor.features = combined.features, verbose = FALSE)
combined_primary_NT.integrated <- IntegrateData(anchorset = SCtransform.anchors, normalization.method = "SCT",
                                          verbose = FALSE)


save(combined_primary_NT.integrated, file=here("data","combined_primary_NT_integrated.RData"))



load(here("data","combined_primary_NT_integrated.RData"))


combined_primary_NT.integrated <- RunPCA(combined_primary_NT.integrated, verbose = FALSE)
combined_primary_NT.integrated <- RunUMAP(combined_primary_NT.integrated, dims = 1:20)
combined_primary_NT.integrated <- RunTSNE(combined_primary_NT.integrated, dims = 1:20)
combined_primary_NT.integrated <- FindNeighbors(combined_primary_NT.integrated, reduction = "pca", dims = 1:20)
combined_primary_NT.integrated <- FindClusters(combined_primary_NT.integrated, resolution = 0.4)


DimPlot(combined_primary_NT.integrated, reduction = "tsne", pt.size=1, label=T)
DimPlot(combined_primary_NT.integrated, reduction = "umap", pt.size=1, label=T)

DimPlot(combined_primary_NT.integrated, reduction = "umap", split.by = "orig.ident", pt.size=1, label=T)
DimPlot(combined_primary_NT.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=1)

load(here("output","primary_predict_transfer.RData"))
combined_primary_NT.integrated <- AddMetaData(combined_primary_NT.integrated, metadata = pri_pred)



combined_primary_NT.integrated@meta.data$cluster_ID[which(is.na(combined_primary_NT.integrated@meta.data$cluster_ID))]<-"organoid"
combined_primary_NT.integrated@meta.data$cluster_ID<-as.factor(combined_primary_NT.integrated@meta.data$cluster_ID)
combined_primary_NT.integrated@meta.data$cluster_ID<-factor(combined_primary_NT.integrated@meta.data$cluster_ID, levels = c(
  "BEST4 enterocyte","crypt","enterocyte",
  "enteroendocrine","Goblet cell", 
  "Paneth","Paneth (UC only)",
  "TA","Tuft","organoid"))
cols_manual_less<-c("#3D0AF2","#6521C1","#367C34",
                     "#B921C1","#308AC8",
                    "#BFC121","#CFA218",
                    "#C86E30","#C12134","lightgrey")


grid.arrange(
  DimPlot(combined_primary_NT.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=cols_manual_less),
  DimPlot(combined_primary_NT.integrated, reduction = "umap", group.by = "integrated_snn_res.0.25", pt.size=0.5, label=T), ncol=2, widths=c(1,0.75)) 


ggsave(file=here("figs","NT_integrated_with_primary.pdf"),
       grid.arrange(
         DimPlot(combined_primary_NT.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=cols_manual_less),
         DimPlot(combined_primary_NT.integrated, reduction = "umap", group.by = "integrated_snn_res.0.25", pt.size=0.5, label=T), ncol=2, widths=c(1,0.75)) ,
       w=11,h=5)
ggsave(file=here("figs/jpeg","NT_integrated_with_primary.jpeg"),
       grid.arrange(
         DimPlot(combined_primary_NT.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=cols_manual_less),
         DimPlot(combined_primary_NT.integrated, reduction = "umap", group.by = "integrated_snn_res.0.25", pt.size=0.5, label=T), ncol=2, widths=c(1,0.75)),
       w=11,h=5)



  
#'## Relabel NT clusters
NT.cells.integrated
NT.cells.integrated$cluster_ID<-"Unidentified"
NT.cells.integrated$cluster_ID[which(NT.cells.integrated$seurat_clusters==0)]<-"enterocyte"
NT.cells.integrated$cluster_ID[which(NT.cells.integrated$seurat_clusters==1)]<-"crypt"
NT.cells.integrated$cluster_ID[which(NT.cells.integrated$seurat_clusters==2)]<-"crypt"
NT.cells.integrated$cluster_ID[which(NT.cells.integrated$seurat_clusters==3)]<-"crypt"
NT.cells.integrated$cluster_ID[which(NT.cells.integrated$seurat_clusters==4)]<-"enterocyte"
NT.cells.integrated$cluster_ID[which(NT.cells.integrated$seurat_clusters==5)]<-"TA"
NT.cells.integrated$cluster_ID[which(NT.cells.integrated$seurat_clusters==6)]<-"TA"


DimPlot(NT.cells.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=c("#6baed6","#238b45","#f16913"))
ggsave(file=here("figs","UMAP_NT_labelled.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","UMAP_NT_labelled.jpeg"), w=6,h=4)

DimPlot(NT.cells.integrated, reduction = "pca", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=c("#6baed6","#238b45","#f16913"))
ggsave(file=here("figs","PCA_NT_labelled.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","PCA_NT_labelled.jpeg"), w=6,h=4)


 


##################################
#' # IFNg
##################################
IFNg.cells.integrated<-readRDS(here("data","IFNg_normalized_integrated.rds"))


###########
#' ## Seurat transfer
###########
org_pri.anchors <- FindTransferAnchors(reference = primary.epi.cells_subset, query = IFNg.cells.integrated, dims = 1:20, project.query = TRUE )
predictions_primary <- TransferData(anchorset = org_pri.anchors, refdata = primary.epi.cells_subset$cluster_ID, dims = 1:20, k.weight=10)
org_pri.cells.integrated <- AddMetaData(IFNg.cells.integrated, metadata = predictions_primary)

#'  Bar plot
cnt_plt<-as.data.frame(table(org_pri.cells.integrated$seurat_clusters, org_pri.cells.integrated$predicted.id))

ggplot(cnt_plt, aes(x=Var1, y=Freq,fill=Var2 )) + geom_bar(position="stack", stat="identity", color="black")+theme_bw()+th+xlab("Cluster")+ylab("Cell Number")+
  scale_fill_manual(values=c("#fb6a4a","#e31a1c","#2171b5","#238443","#74c476","#addd8e","#6baed6","#ffff33","#00441b"), name="Primary Epi Type")
ggsave(file=here("figs","IFNg_seurat_transfer_clusters_bar.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","IFNg_seurat_transfer_clusters_bar.jpeg"),w=6,h=4)


###########
#' ## Integration
###########
combined_primary_IFNg<-list(IFNg.cells.integrated, primary.epi.cells_subset)

combined.features <- SelectIntegrationFeatures(object.list = combined_primary_IFNg, nfeatures = 3000)

options(future.globals.maxSize = 3000 * 1024^2)
combined_primary_IFNg <- PrepSCTIntegration(object.list = combined_primary_IFNg, anchor.features = combined.features,verbose = FALSE)

SCtransform.anchors <- FindIntegrationAnchors(object.list = combined_primary_IFNg, normalization.method = "SCT",
                                              anchor.features = combined.features, verbose = FALSE)
combined_primary_IFNg.integrated <- IntegrateData(anchorset = SCtransform.anchors, normalization.method = "SCT",
                                                verbose = FALSE)


save(combined_primary_IFNg.integrated, file=here("data","combined_primary_IFNg_integrated.RData"))



load(here("data","combined_primary_IFNg_integrated.RData"))


combined_primary_IFNg.integrated <- RunPCA(combined_primary_IFNg.integrated, verbose = FALSE)
combined_primary_IFNg.integrated <- RunUMAP(combined_primary_IFNg.integrated, dims = 1:20)
combined_primary_IFNg.integrated <- RunTSNE(combined_primary_IFNg.integrated, dims = 1:20)
combined_primary_IFNg.integrated <- FindNeighbors(combined_primary_IFNg.integrated, reduction = "pca", dims = 1:20)
combined_primary_IFNg.integrated <- FindClusters(combined_primary_IFNg.integrated, resolution = 0.4)


DimPlot(combined_primary_IFNg.integrated, reduction = "tsne", pt.size=1, label=T)
DimPlot(combined_primary_IFNg.integrated, reduction = "umap", pt.size=1, label=T)

DimPlot(combined_primary_IFNg.integrated, reduction = "umap", split.by = "orig.ident", pt.size=1, label=T)
DimPlot(combined_primary_IFNg.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=1)

load(here("output","primary_predict_transfer.RData"))
combined_primary_IFNg.integrated <- AddMetaData(combined_primary_IFNg.integrated, metadata = pri_pred)



combined_primary_IFNg.integrated@meta.data$cluster_ID[which(is.na(combined_primary_IFNg.integrated@meta.data$cluster_ID))]<-"organoid"
combined_primary_IFNg.integrated@meta.data$cluster_ID<-as.factor(combined_primary_IFNg.integrated@meta.data$cluster_ID)
combined_primary_IFNg.integrated@meta.data$cluster_ID<-factor(combined_primary_IFNg.integrated@meta.data$cluster_ID, levels = c(
  "BEST4 enterocyte","crypt","enterocyte",
  "enteroendocrine","Goblet cell", 
  "Paneth","Paneth (UC only)",
  "TA","Tuft","organoid"))
cols_manual_less<-c("#3D0AF2","#6521C1","#367C34",
                    "#B921C1","#308AC8",
                    "#BFC121","#CFA218",
                    "#C86E30","#C12134","lightgrey")

grid.arrange(
  DimPlot(combined_primary_IFNg.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=cols_manual_less),
  DimPlot(combined_primary_IFNg.integrated, reduction = "umap", group.by = "integrated_snn_res.0.25", pt.size=0.5, label=T), ncol=2, widths=c(1,0.75))



ggsave(file=here("figs","IFNg_integrated_with_primary.pdf"),
       grid.arrange(
         DimPlot(combined_primary_IFNg.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=cols_manual_less),
         DimPlot(combined_primary_IFNg.integrated, reduction = "umap", group.by = "integrated_snn_res.0.25", pt.size=0.5, label=T), ncol=2, widths=c(1,0.75)) ,
       w=11,h=5)
ggsave(file=here("figs/jpeg","IFNg_integrated_with_primary.jpeg"),
       grid.arrange(
         DimPlot(combined_primary_IFNg.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=cols_manual_less),
         DimPlot(combined_primary_IFNg.integrated, reduction = "umap", group.by = "integrated_snn_res.0.25", pt.size=0.5, label=T), ncol=2, widths=c(1,0.75)),
       w=11,h=5)


#'## Relabel IFNg clusters
IFNg.cells.integrated
IFNg.cells.integrated$cluster_ID<-"Unidentified"
IFNg.cells.integrated$cluster_ID[which(IFNg.cells.integrated$seurat_clusters==0)]<-"TA"
IFNg.cells.integrated$cluster_ID[which(IFNg.cells.integrated$seurat_clusters==1)]<-"crypt"
IFNg.cells.integrated$cluster_ID[which(IFNg.cells.integrated$seurat_clusters==2)]<-"enterocyte"
IFNg.cells.integrated$cluster_ID[which(IFNg.cells.integrated$seurat_clusters==3)]<-"crypt"
IFNg.cells.integrated$cluster_ID[which(IFNg.cells.integrated$seurat_clusters==4)]<-"enterocyte"
IFNg.cells.integrated$cluster_ID[which(IFNg.cells.integrated$seurat_clusters==5)]<-"TA"


DimPlot(IFNg.cells.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=c("#6baed6","#238b45","#f16913","pink"))
ggsave(file=here("figs","UMAP_IFNg_labelled.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","UMAP_IFNg_labelled.jpeg"), w=6,h=4)

DimPlot(IFNg.cells.integrated, reduction = "pca", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=c("#6baed6","#238b45","#f16913","pink"))
ggsave(file=here("figs","PCA_IFNg_labelled.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","PCA_IFNg_labelled.jpeg"), w=6,h=4)





########################################
#' # TNFa organoids
########################################
TNFa.cells.integrated<-readRDS(here("data","TNFa_normalized_integrated.rds"))




###########
#' ## Seurat transfer
###########
org_pri.anchors <- FindTransferAnchors(reference = primary.epi.cells_subset, query = TNFa.cells.integrated, dims = 1:20, project.query = TRUE )
predictions_primary <- TransferData(anchorset = org_pri.anchors, refdata = primary.epi.cells_subset$cluster_ID, dims = 1:20, k.weight=10)
org_pri.cells.integrated <- AddMetaData(TNFa.cells.integrated, metadata = predictions_primary)

#'  Bar plot
cnt_plt<-as.data.frame(table(org_pri.cells.integrated$seurat_clusters, org_pri.cells.integrated$predicted.id))

ggplot(cnt_plt, aes(x=Var1, y=Freq,fill=Var2 )) + geom_bar(position="stack", stat="identity", color="black")+theme_bw()+th+xlab("Cluster")+ylab("Cell Number")+
  scale_fill_manual(values=c("#fb6a4a","#e31a1c","#2171b5","#238443","#74c476","#addd8e","#6baed6","#ffff33","#00441b"), name="Primary Epi Type")
ggsave(file=here("figs","TNFa_seurat_transfer_clusters_bar.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","TNFa_seurat_transfer_clusters_bar.jpeg"),w=6,h=4)

###########
#' ## Integration together
###########
combined_primary_TNFa<-list(TNFa.cells.integrated, primary.epi.cells_subset)

combined.features <- SelectIntegrationFeatures(object.list = combined_primary_TNFa, nfeatures = 3000)

options(future.globals.maxSize = 3000 * 1024^2)
combined_primary_TNFa <- PrepSCTIntegration(object.list = combined_primary_TNFa, anchor.features = combined.features,verbose = FALSE)

SCtransform.anchors <- FindIntegrationAnchors(object.list = combined_primary_TNFa, normalization.method = "SCT",
                                              anchor.features = combined.features, verbose = FALSE)
combined_primary_TNFa.integrated <- IntegrateData(anchorset = SCtransform.anchors, normalization.method = "SCT",
                                                verbose = FALSE)


save(combined_primary_TNFa.integrated, file=here("data","combined_primary_TNFa_integrated.RData"))



load(here("data","combined_primary_TNFa_integrated.RData"))


combined_primary_TNFa.integrated <- RunPCA(combined_primary_TNFa.integrated, verbose = FALSE)
combined_primary_TNFa.integrated <- RunUMAP(combined_primary_TNFa.integrated, dims = 1:20)
combined_primary_TNFa.integrated <- RunTSNE(combined_primary_TNFa.integrated, dims = 1:20)
combined_primary_TNFa.integrated <- FindNeighbors(combined_primary_TNFa.integrated, reduction = "pca", dims = 1:20)
combined_primary_TNFa.integrated <- FindClusters(combined_primary_TNFa.integrated, resolution = 0.4)


DimPlot(combined_primary_TNFa.integrated, reduction = "tsne", pt.size=1, label=T)
DimPlot(combined_primary_TNFa.integrated, reduction = "umap", pt.size=1, label=T)

DimPlot(combined_primary_TNFa.integrated, reduction = "umap", split.by = "orig.ident", pt.size=1, label=T)
DimPlot(combined_primary_TNFa.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=1)

load(here("output","primary_predict_transfer.RData"))
combined_primary_TNFa.integrated <- AddMetaData(combined_primary_TNFa.integrated, metadata = pri_pred)



combined_primary_TNFa.integrated@meta.data$cluster_ID[which(is.na(combined_primary_TNFa.integrated@meta.data$cluster_ID))]<-"organoid"
combined_primary_TNFa.integrated@meta.data$cluster_ID<-as.factor(combined_primary_TNFa.integrated@meta.data$cluster_ID)
combined_primary_TNFa.integrated@meta.data$cluster_ID<-factor(combined_primary_TNFa.integrated@meta.data$cluster_ID, levels = c(
  "BEST4 enterocyte","crypt","enterocyte",
  "enteroendocrine","Goblet cell", 
  "Paneth","Paneth (UC only)",
  "TA","Tuft","organoid"))
cols_manual_less<-c("#3D0AF2","#6521C1","#367C34",
                    "#B921C1","#308AC8",
                    "#BFC121","#CFA218",
                    "#C86E30","#C12134","lightgrey")



grid.arrange(
  DimPlot(combined_primary_TNFa.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=cols_manual_less),
  DimPlot(combined_primary_TNFa.integrated, reduction = "umap", group.by = "integrated_snn_res.0.25", pt.size=0.5, label=T), ncol=2, widths=c(1,0.75))


# 
ggsave(file=here("figs","TNFa_integrated_with_primary.pdf"),
       grid.arrange(
         DimPlot(combined_primary_TNFa.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=cols_manual_less),
         DimPlot(combined_primary_TNFa.integrated, reduction = "umap", group.by = "integrated_snn_res.0.25", pt.size=0.5, label=T), ncol=2, widths=c(1,0.75)) ,
       w=11,h=5)
ggsave(file=here("figs/jpeg","TNFa_integrated_with_primary.jpeg"),
       grid.arrange(
         DimPlot(combined_primary_TNFa.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=cols_manual_less),
         DimPlot(combined_primary_TNFa.integrated, reduction = "umap", group.by = "integrated_snn_res.0.25", pt.size=0.5, label=T), ncol=2, widths=c(1,0.75)),
       w=11,h=5)


#'## Relabel TNFa clusters
TNFa.cells.integrated
TNFa.cells.integrated$cluster_ID<-"Unidentified"
TNFa.cells.integrated$cluster_ID[which(TNFa.cells.integrated$seurat_clusters==0)]<-"enterocyte"
TNFa.cells.integrated$cluster_ID[which(TNFa.cells.integrated$seurat_clusters==1)]<-"crypt"
TNFa.cells.integrated$cluster_ID[which(TNFa.cells.integrated$seurat_clusters==2)]<-"enterocyte"
TNFa.cells.integrated$cluster_ID[which(TNFa.cells.integrated$seurat_clusters==3)]<-"crypt"
TNFa.cells.integrated$cluster_ID[which(TNFa.cells.integrated$seurat_clusters==4)]<-"crypt"
TNFa.cells.integrated$cluster_ID[which(TNFa.cells.integrated$seurat_clusters==5)]<-"TA"


DimPlot(TNFa.cells.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=c("#6baed6","#238b45","#f16913","pink"))
ggsave(file=here("figs","UMAP_TNFa_labelled.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","UMAP_TNFa_labelled.jpeg"), w=6,h=4)

DimPlot(TNFa.cells.integrated, reduction = "pca", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=c("#6baed6","#238b45","#f16913","pink"))
ggsave(file=here("figs","PCA_TNFa_labelled.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","PCA_TNFa_labelled.jpeg"), w=6,h=4)


save(TNFa.cells.integrated, IFNg.cells.integrated, NT.cells.integrated, file=here("data","organoids_labelled_stimulations.RData"))


TNF_celllabel<-TNFa.cells.integrated@meta.data[,c("orig.ident","cluster_ID")]
IFNg_celllabel<-IFNg.cells.integrated@meta.data[,c("orig.ident","cluster_ID")]
NT_celllabel<-NT.cells.integrated@meta.data[,c("orig.ident","cluster_ID")]

cell_typelabel<-rbind(NT_celllabel, IFNg_celllabel,TNF_celllabel)
save(cell_typelabel, file=here("output","celltype_labels_organoid.RData"))



#' 
#' 
#' ###########
#' #' ## Cell proportions
#' ###########
#' d10x.stim<-readRDS(here("data","d10x_raw_merged.rds"))
#' 
#' ## add cell type labels from split analysis
#' load(here("output","celltype_labels_organoid.RData"))
#' identical(rownames(cell_typelabel),rownames(d10x.stim@meta.data))
#' cell_typelabel$cluster_ID[which(cell_typelabel$cluster_ID=="entero_stem")]<-"enterocyte_stem"
#' d10x.stim <- AddMetaData(d10x.stim, metadata = cell_typelabel)
#' 
#' 
#' proportions<-do.call(rbind, lapply(unique(d10x.stim@meta.data$individual), function(x){
#'   donor<-d10x.stim@meta.data[which(d10x.stim@meta.data$individual==x),]
#'   count<-as.data.frame(tapply(rownames(donor), list(donor$orig.ident, donor$cluster_ID), length))
#'   count<-count/rowSums(count)
#'   count$stimulation<-rownames(count)
#'   count$donor<-x
#'   count
#' }))
#' 
#' 
#' 
#' proportions
#' 
#' 
#' plt_pro<-melt(proportions)
#' plt_pro$stimulation<-factor(plt_pro$stimulation, levels=c("NT","IFNg","TNFa"))
#' 
#' ggplot(plt_pro, aes(stimulation, value, fill=stimulation))+geom_boxplot(outlier.shape=NA)+facet_wrap(~variable)+geom_line(aes(stimulation, value, group=donor), color="grey")+
#'   geom_point(shape=21, color="black", size=2)+theme_bw()+th+scale_fill_manual(values=c("grey80","cornflowerblue","firebrick4"))
#' ggsave(file=here("figs","proportion_celltype.pdf"), w=6,h=4)
#' ggsave(file=here("figs/jpeg","proportion_celltype.jpeg"), w=6,h=4)
#' 
#' 
#' 
#' ## univariate paired test in each cell type comparing NT to stimulations
#'   do.call(rbind, lapply(c("crypt","enterocyte","TA"), function(cell){
#'   data.frame(cell=cell,
#'              IFNg_pval=wilcox.test(proportions[which(proportions$stimulation=="NT"),cell], proportions[which(proportions$stimulation=="IFNg"),cell], paired = TRUE, alternative = "two.sided")$p.value,
#'              TNFa_pval=wilcox.test(proportions[which(proportions$stimulation=="NT"),cell], proportions[which(proportions$stimulation=="TNFa"),cell], paired = TRUE, alternative = "two.sided")$p.value)}))
#' 
#' ## univariate unpaired test in each cell type comparing NT to stimulations
#'  wil_p<- do.call(rbind,lapply(c("crypt","enterocyte","TA"), function(cell){
#'     data.frame(cell=cell,
#'                IFNg_pval=wilcox.test(proportions[which(proportions$stimulation=="NT"),cell], proportions[which(proportions$stimulation=="IFNg"),cell], paired = F, alternative = "two.sided")$p.value,
#'                TNFa_pval=wilcox.test(proportions[which(proportions$stimulation=="NT"),cell], proportions[which(proportions$stimulation=="TNFa"),cell], paired = F, alternative = "two.sided")$p.value)}))
#'   
#'   p.adjust(wil_p[,2], method="fdr", n=6)
#'   p.adjust(wil_p[,3], method="fdr", n=6)
#'   


#'## R Session Info
sessionInfo()
  