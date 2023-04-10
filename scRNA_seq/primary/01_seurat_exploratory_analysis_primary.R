
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

options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))



##############################################
#' PRIMARY
##############################################
dataset_loc <- here("../ibd/data/raw/scRNAseq/scRNASEQ_DATA/Pediatric_samples/pediatric_cellranger302_GRCh38-3.0.0")

meta_primary<-read.csv(here("../ibd/data/raw/scRNAseq/scRNASEQ_DATA/Pediatric_samples/Pediatric_sample_tracker.csv"))
colnames(meta_primary)[2]<-"donor"

## From a TRIPP sheet Rasa shared
meta_primary_update<-read.csv(here("data","primary_sc_TRIPP_update.csv"))
identical(meta_primary_update$Sample.name, meta_primary$donor)
meta_primary$Diagnosis<-c(meta_primary_update$Diagnosis_update)
meta_primary$Diagnosis[which(is.na(meta_primary$Diagnosis))]<-"neo"


#' Cohort table
diagnosis<-c("CD","Ctrl","UC","neo")

d10x.list_primary<-lapply(1:length(diagnosis),function(y) {
  ids <- as.character(meta_primary$Sanger.Sample.ID)[grep(diagnosis[y],meta_primary$Diagnosis)]
  
  d10x.data <- sapply(ids, function(i){
    d10x <- Read10X(file.path(dataset_loc,paste("outs_",i, sep=""),"filtered_feature_bc_matrix"))
    colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
    d10x
  })
  
  d10x.data<-do.call("cbind",d10x.data)
  
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x <- CreateSeuratObject(counts = d10x.data, project = diagnosis[y], min.cells = 3, min.features = 0)
  
  meta_cell<-data.frame(cell=colnames(d10x), individual=sapply(colnames(d10x), function(x) strsplit(x,"-")[[1]][2]))
  meta_cell_add<-merge(meta_cell, meta_primary, by.x="individual", by.y="Sanger.Sample.ID")
  #meta_cell_add<-meta_cell_add[match(meta_cell_add$cell, colnames(d10x)),]
  meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
  
  identical(meta_cell_add$cell, colnames(d10x))
  
  d10x<- AddMetaData(d10x, meta_cell_add$donor, col.name = "individual")})

names(d10x.list_primary)<-diagnosis


## cell counts
plt_count_raw<-do.call(rbind,lapply(1:length(d10x.list_primary), function(x) {
  df<-as.data.frame(tapply(rownames(d10x.list_primary[[x]]@meta.data), list(d10x.list_primary[[x]]@meta.data$individual, d10x.list_primary[[x]]@meta.data$orig.ident), length))
  df$individual<-rownames(df)
  df$condition<-names(d10x.list_primary)[x]
  colnames(df)<-c("cell_count","individual","condition")
  df}))
plt_count_raw

ggplot(plt_count_raw, aes(condition, cell_count))+
  geom_boxplot(fill="lightgrey")+geom_point()+
  theme_bw()+geom_text(aes(label=individual), hjust=-0.25)+ylim(0,15000)+ylab("Total Cell Number")



#'## QC
#'The percentage of reads that map to the mitochondrial genome
#'Low-quality / dying cells often exhibit extensive mitochondrial contamination
#'We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
#'We use the set of all genes starting with MT- as a set of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
invisible(lapply(1:length(d10x.list_primary), function(x){
  d10x.list_primary[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list_primary[[x]], pattern = "^MT-")}))

# Show QC metrics for the first 5 cells
head(d10x.list_primary[[2]]@meta.data, 5)

#'Low-quality cells or empty droplets will often have very few genes
#'Cell doublets or multiplets may exhibit an aberrantly high gene count


# Visualize QC metrics
#nFeature number of unique genes
#nCount number of total molecules
plt_QC_data<-do.call(rbind, lapply(1:4, function(x) d10x.list_primary[[x]]@meta.data))

ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

ggsave(file=here("figs","primary_ncount_nfeatures.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","primary_ncount_nfeatures.jpeg"), w=6,h=4)


ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 10)+ theme_bw()+xlab("Percent Mitochondrial")+th

ggsave(file=here("figs","primary_MT_features.pdf"), w=4,h=2)
ggsave(file=here("figs/jpeg","primary_MT_features.jpeg"), w=4,h=2)




#' #### scrublet results
samples<-meta_primary$Sanger.Sample.ID

doublet_scores_all<-lapply(1:length(samples), function(x){
  ## scrublet output
  doublet_scores<-read.csv(paste(here("data/"),"scrublet_output_table_",samples[x],".csv",sep=""))
  doublet_scores$index_edit<-gsub("1",samples[x],doublet_scores$index)
  doublet_scores
})
doublet_scores_all<-do.call(rbind, doublet_scores_all)

## Match doublet information
doublet_scores_list<-lapply(1:length(d10x.list_primary), function(x){
  doublet_scores_all<-doublet_scores_all[match(colnames(d10x.list_primary[[x]]), doublet_scores_all$index_edit),]
  identical(colnames(d10x.list_primary[[x]]), doublet_scores_all$index_edit)
  doublet_scores_all
})

invisible(lapply(1:length(d10x.list_primary), function(x){
  d10x.list_primary[[x]]<<- AddMetaData(d10x.list_primary[[x]], doublet_scores_list[[x]]$doublet_score, col.name = "doublet_score")
  d10x.list_primary[[x]]<<- AddMetaData(d10x.list_primary[[x]], doublet_scores_list[[x]]$predicted_doublet, col.name = "predicted_doublet")
}))

plt_QC_data<-do.call(rbind, lapply(1:4, function(x) d10x.list_primary[[x]]@meta.data))

ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=doublet_score)) + 
  geom_point(size=0.75) + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Doublet\nScore") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

ggsave(file=here("figs","primary_doublet_features.pdf"),w=5,h=3)
ggsave(file=here("figs/jpeg","primary_doublet_features.jpeg"), w=5,h=3)

ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=doublet_score)) + 
  geom_point(size=0.75) + facet_wrap(~predicted_doublet)+
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Doublet\nScore") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

ggsave(file=here("figs","primary_doublet_features_split.pdf"), w=12,h=4)
ggsave(file=here("figs/jpeg","primary_doublet_features_spit.jpeg"), w=12,h=4)

table(plt_QC_data$predicted_doublet)
dim(plt_QC_data)
table(plt_QC_data$predicted_doublet, plt_QC_data$individual)



## what is 3 MAD from mean
median(plt_QC_data$nFeature_RNA)+(mad(plt_QC_data$nFeature_RNA)*3)
median(plt_QC_data$nFeature_RNA)-(mad(plt_QC_data$nFeature_RNA)*3)


#'We filter cells that have unique feature counts over 6,000 or less than 500
#'We filter cells that have >10% mitochondrial counts
#'we will also filter doublets as called by scrublet
d10x.list_primary.raw<-d10x.list_primary

invisible(lapply(1:length(d10x.list_primary), function(x){
  d10x.list_primary[[x]] <<- subset(d10x.list_primary[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 & predicted_doublet=="False")
}))

d10x.list_primary

## cell counts after QC
plt_count_QC<-do.call(rbind,lapply(1:length(d10x.list_primary), function(x) {
  df<-as.data.frame(tapply(rownames(d10x.list_primary[[x]]@meta.data), list(d10x.list_primary[[x]]@meta.data$individual, d10x.list_primary[[x]]@meta.data$orig.ident), length))
  df$individual<-rownames(df)
  df$condition<-names(d10x.list_primary)[x]
  colnames(df)<-c("cell_count","individual","condition")
  df}))
plt_count_QC

#save to plot fancy
save(plt_count_raw, plt_count_QC, file=here("data","primary_cell_count.RData"))

plt_count_raw<-plt_count_raw[which(plt_count_raw$condition!="neo"),]
plt_count_raw$condition<-as.factor(plt_count_raw$condition)
levels(plt_count_raw$condition)<-c("CD","Control","UC")
plt_count_raw$condition<-factor(plt_count_raw$condition, levels=c("Control","CD","UC"))

plt_count_QC<-plt_count_QC[which(plt_count_QC$condition!="neo"),]
plt_count_QC$condition<-as.factor(plt_count_QC$condition)
levels(plt_count_QC$condition)<-c("CD","Control","UC")
plt_count_QC$condition<-factor(plt_count_QC$condition, levels=c("Control","CD","UC"))


## only samples in thesis
plt_count_raw<-plt_count_raw[-grep("036NEG|036POS",plt_count_raw$individual),]
plt_count_QC<-plt_count_QC[-grep("036NEG|036POS",plt_count_QC$individual),]

cell_count<-grid.arrange(ggplot(plt_count_raw, aes(condition, cell_count,fill=condition))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=individual), hjust=-0.25, size=3)+ylim(0, 15000)+xlab("Diagnosis")+
                           ylab("Total Cell Number")+th+fillscale_diagnosis+
                           theme(legend.position = "none")+ggtitle("Before Quality Control"),
                         ggplot(plt_count_QC, aes(condition, cell_count,fill=condition))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=individual), hjust=-0.25, size=3)+ylim(0, 15000)+xlab("Diagnosis")+
                           ylab("Total Cell Number")+th+fillscale_diagnosis+
                           theme(legend.position = "none")+ggtitle("After Quality Control"), ncol=2)

ggsave(cell_count,file=here("figs","primary_scCount.pdf"), w=8,h=4)
ggsave(cell_count,file=here("figs/jpeg","primary_scCount.jpeg"), w=8,h=4)


# $CD
# An object of class Seurat 
# 23944 features across 24996 samples within 1 assay 
# Active assay: RNA (23944 features, 0 variable features)
# 
# $Ctrl
# An object of class Seurat 
# 23904 features across 27622 samples within 1 assay 
# Active assay: RNA (23904 features, 0 variable features)
# 
# $UC
# An object of class Seurat 
# 19768 features across 2703 samples within 1 assay 
# Active assay: RNA (19768 features, 0 variable features)
# 
# $neo
# An object of class Seurat 
# 20196 features across 3743 samples within 1 assay 
# Active assay: RNA (20196 features, 0 variable features)


d10x.primary <- merge(d10x.list_primary[[1]], d10x.list_primary[[2]])
d10x.primary <- merge(d10x.primary, d10x.list_primary[[3]])
d10x.primary <- merge(d10x.primary, d10x.list_primary[[4]])

              # An object of class Seurat 
              # 24939 features across 59320 samples within 1 assay 
              # Active assay: RNA (24939 features, 0 variable features)

d10x.primary.raw<-d10x.primary

saveRDS(d10x.primary.raw, file = here("data","d10x_primary_raw_merged.rds"))



d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))



d10x.primary <- NormalizeData(d10x.primary)
d10x.primary <- FindVariableFeatures(d10x.primary, selection.method = "vst", nfeatures = 2000)
d10x.primary <- ScaleData(d10x.primary) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
d10x.primary <- RunPCA(d10x.primary, ndims.print = 1:10, nfeatures.print = 10)
d10x.primary <- RunUMAP(d10x.primary, dims = 1:30)
d10x.primary <- RunTSNE(d10x.primary, dims = 1:30)



######################
## cell cycle gene expression
######################

## which PC associated to nfeature? 2
pca_mat<-as.data.frame(Embeddings(object = d10x.primary, reduction = "pca"))
pca_mat$cell<-rownames(pca_mat)
meta_primary<-d10x.primary@meta.data
meta_primary$cell<-rownames(meta_primary)
plt<-merge(pca_mat, meta_primary,by="cell")

sapply(2:51, function(x) cor(plt[,x],plt[,"nFeature_RNA"]))
FeaturePlot(d10x.primary, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
## PC2 tops genes include cell cycle genes: "TPX2"   "TOP2A"  "MKI67"  "CENPF"  "SMC4"   "TUBB4B"

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x.primary <- CellCycleScoring(d10x.primary, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

DimPlot(d10x.primary, reduction="pca")
ggsave(file=here("figs","primary_cell_cycle_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","primary_cell_cycle_PCA.jpeg"), w=5,h=5)

FeaturePlot(d10x.primary, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","primary_nfeature_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","primary_nfeature_PCA.jpeg"), w=5.25,h=5)

DimPlot(d10x.primary, reduction="umap", group.by = "Phase")


## regress out cell cycle and other covariates
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#By default, sctransform accounts for cellular sequencing depth, or nUMIs.
d10x.primary <- SCTransform(d10x.primary, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE)


# dimension reduction
d10x.primary <- RunPCA(d10x.primary, verbose = FALSE)
d10x.primary <- RunUMAP(d10x.primary, dims = 1:30)
d10x.primary <- RunTSNE(d10x.primary, dims = 1:30)

# cluster
d10x.primary <- FindNeighbors(d10x.primary, reduction = "pca", dims = 1:20)
d10x.primary <- FindClusters(d10x.primary, resolution = 0.6)

# relabel meta data
d10x.primary@meta.data$orig.ident[which(d10x.primary@meta.data$orig.ident=="Ctrl")]<-"Control"
d10x.primary@meta.data$orig.ident[which(d10x.primary@meta.data$orig.ident=="neo")]<-"Neonatal"

saveRDS(d10x.primary, file = here("data","d10x_primary_normalized.rds"))





d10x.primary<-readRDS(here("data","d10x_primary_normalized.rds"))
# An object of class Seurat 
# 48071 features across 59064 samples within 2 assays 
# Active assay: SCT (23126 features, 3000 variable features)
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, umap, tsne

## visualize
DimPlot(d10x.primary, reduction = "umap", group.by = "Phase", pt.size=1)
DimPlot(d10x.primary, reduction = "umap", group.by = "Phase",split.by = "orig.ident", pt.size=1)

DimPlot(d10x.primary, reduction = "umap", pt.size=0.25, label=T)
ggsave(file=here("figs","primary_UMAP_cluster.pdf"), w=6,h=5)
ggsave(file=here("figs/jpeg","primary_UMAP_cluster.jpeg"), w=6,h=5)

DimPlot(d10x.primary, reduction = "tsne", pt.size=0.25, label=T)
ggsave(file=here("figs","primary_tsne_cluster.pdf"), w=6,h=5)
ggsave(file=here("figs/jpeg","primary_tsne_cluster.jpeg"), w=6,h=5)


DimPlot(d10x.primary, reduction = "umap", group.by = "orig.ident", pt.size=0.25)+colscale_diagnosis
ggsave(file=here("figs","primary_UMAP_raw_stimulation.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","primary_UMAP_raw_stimulation.jpeg"), w=7,h=6)

DimPlot(d10x.primary, reduction = "tsne", group.by = "orig.ident", pt.size=0.25)+colscale_diagnosis
ggsave(file=here("figs","primary_TSNE_raw_stimulation.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","primary_TSNE_raw_stimulation.jpeg"), w=7,h=6)

DimPlot(d10x.primary, reduction="pca", pt.size=1)

DimPlot(d10x.primary, reduction = "umap", group.by = "individual", pt.size=1)
DimPlot(d10x.primary, reduction = "umap", split.by = "individual", group.by = "orig.ident", pt.size=1)


## Low quality cluster?
FeaturePlot(d10x.primary, features = "nCount_RNA", min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.primary, features = "nCount_RNA",split.by = "orig.ident", min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.primary, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)

FeaturePlot(d10x.primary, features = "doublet_score", pt.size=1)
ggsave(file=here("figs","primary_doublet socre.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","primary_doublet socre.jpeg"), w=5,h=5)

DimPlot(d10x.primary, reduction="pca", group.by = "orig.ident", pt.size=1)

DimPlot(d10x.primary, reduction="pca", group.by="Phase")
ggsave(file=here("figs","primary_cell_cycle_PCA_after_regression.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","primary_cell_cycle_PCA_after_regression.jpeg"), w=5,h=5)

FeaturePlot(d10x.primary, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","primary_nfeature_PCA_after_regression.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","primary_nfeature_PCA_after_regression.jpeg"), w=5.25,h=5)




#########################
## Marker genes
#########################
#log-normalized versions of corrected counts
d10x.primary.exp<-as.data.frame(d10x.primary[["SCT"]]@data)

genes<-unique(c("EPCAM", "KRT8", "KRT18", #epi
                "COL1A1", "COL1A2", "COL6A1", "COL6A2", "VWF", "PLVAP", "CDH5", "S100B", # stromal
                "CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G")) #immune

d10x.primary.exp["EPCAM",1:10]

d10x.primary.exp.GOI<-d10x.primary.exp[genes,]
d10x.primary.exp.GOI$gene<-rownames(d10x.primary.exp.GOI)
d10x.primary.exp.GOI<-melt(d10x.primary.exp.GOI)#

umap_mat<-as.data.frame(Embeddings(object = d10x.primary, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta_primary<-d10x.primary@meta.data
meta_primary$cell<-rownames(meta_primary)

plt<-merge(d10x.primary.exp.GOI, meta_primary,by.x="variable", by.y="cell")
plt<-merge(plt, umap_mat,by.x="variable", by.y="cell")

save(plt, file=here("output","GOI_primary.RData"))



load(file=here("output","GOI_primary.RData"))

gene_plot<-function(gene_name){
  plt_gene<-plt[which(plt$gene==gene_name),]
  plt_gene<-plt_gene[order(plt_gene$value),]
  ggplot(plt_gene, aes(UMAP_1,UMAP_2, color=value))+
    geom_point(size=1.5)+theme_bw()+
    scale_colour_gradient2( low = "#2b8cbe",#blue
                            mid = "#f0f0f0", # grey
                            high = "#eb8423", #oragne
                            midpoint = max(plt_gene[,"value"])/2,
                            na.value = "#226d94",
                            name="Count")+theme_classic()+
    ggtitle(gene_name)+th+theme(plot.title = element_text(size = 16, face = "bold"))
}

#epithelial
gene_plot("EPCAM")
gene_plot("KRT18")
gene_plot("KRT8")
#stromal
gene_plot("COL1A1")
gene_plot("COL1A2")
gene_plot("COL6A1")
gene_plot("COL6A2")

gene_plot("VWF")
gene_plot("PLVAP")
gene_plot("CDH5")
gene_plot("S100B")

#immune     
gene_plot("CD52")
gene_plot("CD2")
gene_plot("CD3D")
gene_plot("CD3G")
gene_plot("CD3E")
gene_plot("CD79A")
gene_plot("CD79B")
gene_plot("CD68")
gene_plot("CD14")
gene_plot("CD16")
gene_plot("CD68")
gene_plot("CD83")
gene_plot("CSF1R")
gene_plot("FCER1G")




grid.arrange(gene_plot("EPCAM"),gene_plot("KRT18"),gene_plot("KRT8"), ncol=1)
grid.arrange(gene_plot("COL1A1"),gene_plot("CDH5"),gene_plot("S100B"), ncol=1)
grid.arrange(gene_plot("CD52"),gene_plot("CD79A"),gene_plot("CD83"), ncol=1)


ggsave(grid.arrange(gene_plot("EPCAM"),gene_plot("KRT18"),gene_plot("KRT8"), ncol=1),
       file=here("figs/jpeg","Primary_epicellmarker.jpg"), w=5,h=12)

ggsave(grid.arrange(gene_plot("COL1A1"),gene_plot("CDH5"),gene_plot("S100B"), ncol=1),
       file=here("figs/jpeg","Primary_stromalcellmarker.jpg"), w=5,h=12)

ggsave(grid.arrange(gene_plot("CD52"),gene_plot("CD79A"),gene_plot("CD83"), ncol=1),
       file=here("figs/jpeg","Primary_immunecellmarker.jpg"), w=5,h=12)

#####################
### cell compartment signature score from smilie
#####################
#Signature scores were calculated as the mean log 2 (TP10K+1) across all genes in the signature. Each cluster was assigned to the compartment of its maximal score

epigenes<-c("EPCAM", "KRT8", "KRT18") #epi
stromgenes<-c("COL1A1", "COL1A2", "COL6A1", "COL6A2", "VWF", "PLVAP", "CDH5", "S100B") # stromal
immunegenes<-c("CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G") #immune

#####
## call individual cells into compartment
#####
plt$variable<-as.character(plt$variable)

plt_epi<-plt[which(plt$gene%in%epigenes),]
epi<-as.data.frame(tapply(plt_epi$value, plt_epi$variable, mean))

plt_strom<-plt[which(plt$gene%in%stromgenes),]
strom<-as.data.frame(tapply(plt_strom$value, plt_strom$variable, mean))

plt_immune<-plt[which(plt$gene%in%immunegenes),]
immune<-as.data.frame(tapply(plt_immune$value, plt_immune$variable, mean))

cell_compartment<-cbind(epi, strom, immune)
colnames(cell_compartment)<-c("epi","strom","immune")
cell_compartment$variable<-rownames(cell_compartment)

plt_summary<-merge(plt[,c(1,4:ncol(plt))][!duplicated(plt[,c(1,4:ncol(plt))]),], cell_compartment, by="variable")


plt_summary$compartment<-sapply(1:nrow(plt_summary), function(x) {
  compart<-c("epi","strom","immune")[which(plt_summary[x,c("epi","strom","immune")] == max(plt_summary[x,c("epi","strom","immune")]))]
  if(length(compart)==1){compart}else{"Unclear"}
  })



ggplot(plt_summary, aes(UMAP_1,UMAP_2, color=compartment))+
  geom_point(size=1.5)+theme_classic()+th+scale_color_manual(values=c("red","blue","green","grey"))




#######
## smilie did the calling at the custer not individual cell level
#######
plt$variable<-as.character(plt$variable)

plt_epi<-plt[which(plt$gene%in%epigenes),]
epi<-as.data.frame(tapply(plt_epi$value, plt_epi$seurat_clusters, mean))

plt_strom<-plt[which(plt$gene%in%stromgenes),]
strom<-as.data.frame(tapply(plt_strom$value, plt_strom$seurat_clusters, mean))

plt_immune<-plt[which(plt$gene%in%immunegenes),]
immune<-as.data.frame(tapply(plt_immune$value, plt_immune$seurat_clusters, mean))

cell_compartment<-cbind(epi, strom, immune)
colnames(cell_compartment)<-c("epi","strom","immune")
cell_compartment$seurat_clusters<-rownames(cell_compartment)

cell_compartment$compartment<-sapply(1:nrow(cell_compartment), function(x) {
  compart<-c("epi","strom","immune")[which(cell_compartment[x,c("epi","strom","immune")] == max(cell_compartment[x,c("epi","strom","immune")]))]
  if(length(compart)==1){compart}else{"Unclear"}
})


plt_summary<-merge(plt[,c(1,4:ncol(plt))][!duplicated(plt[,c(1,4:18)]),], cell_compartment[,c("seurat_clusters","compartment")], by="seurat_clusters")



ggplot(plt_summary, aes(UMAP_1,UMAP_2, color=compartment))+
  geom_point(size=1.5)+theme_classic()+th+scale_color_manual(values=c("red","blue","green","grey"))


cell_compartment[which(cell_compartment$seurat_clusters==24),]
unique(plt_summary[which(plt_summary$compartment=="epi"),"seurat_clusters"])

save(plt_summary, file=here("output","compartment_primary_assignments.RData"))

# unique(plt_summary$orig.ident)
# plt_summart_proportion<-plt_summary[which(plt_summary$orig.ident!="Neonatal"),]
# plt_summart_proportion<-plt_summart_proportion[-grep("NEG|POS", plt_summart_proportion$individual),]


#### dimension reduction by compartment

plt_summary<-plt_summary[match(colnames(d10x.primary), plt_summary$variable),]
identical(colnames(d10x.primary), plt_summary$variable)

d10x.primary<- AddMetaData(d10x.primary, plt_summary$compartment, col.name = "compartment")



DimPlot(d10x.primary, reduction = "umap", group.by = "compartment", pt.size=0.25)+
  scale_color_manual(values=c("#b2182b","#4393c3","#5aae61"))
ggsave(file=here("figs","primary_UMAP_raw_compartment.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","primary_UMAP_raw_compartment.jpeg"), w=7,h=6)

DimPlot(d10x.primary, reduction = "tsne", group.by = "compartment", pt.size=0.25)+
  scale_color_manual(values=c("#b2182b","#4393c3","#5aae61"))
ggsave(file=here("figs","primary_TSNE_raw_compartment.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","primary_TSNE_raw_compartment.jpeg"), w=7,h=6)

DimPlot(d10x.primary, reduction = "tsne", label=T, pt.size=0.25)
ggsave(file=here("figs","primary_TSNE_raw_cluster.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","primary_TSNE_raw_cluster.jpeg"), w=7,h=6)

DimPlot(d10x.primary, reduction = "umap", label=T, pt.size=0.25)
ggsave(file=here("figs","primary_UMAP_raw_cluster.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","primary_UMAP_raw_cluster.jpeg"), w=7,h=6)



## proportion compartment by diagnosis
donor_counts<-do.call(rbind,lapply(1:length(unique(plt_summary$individual)), function(x){
  donor_count<-plt_summary[which(plt_summary$individual==unique(plt_summary$individual)[x]),]
  counts<-as.data.frame(t(as.data.frame((tapply(donor_count$variable, donor_count$compartment, length)/nrow(donor_count))*100)))
  counts$individual<-unique(plt_summary$individual)[x]
  counts$diagnosis<-unique(plt_summary[which(plt_summary$individual==unique(plt_summary$individual)[x]),"orig.ident"])
   counts
}))
rownames(donor_counts)<-NULL

donor_counts<-melt(donor_counts)

ggplot() + geom_bar(aes(y = value, x = individual, fill = variable), data = donor_counts, color="black",stat="identity")+
  facet_grid(variable~diagnosis, scales="free_x", space="free_x")+theme_bw()+th+
  scale_fill_manual(values=c("#b2182b","#4393c3","#5aae61"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(file=here("figs","primary_composistion.pdf"), w=10,h=4)
ggsave(file=here("figs/jpeg","primary_composistion.jpeg"), w=10,h=4)


ggplot(donor_counts, aes(diagnosis, value, fill=diagnosis))+geom_boxplot(outlier.shape=NA)+
  geom_point(shape=21,size=2,position=position_jitter(w=0.1))+facet_wrap(~variable)+theme_bw()+th+fillscale_diagnosis
ggsave(file=here("figs","primary_composistion_box.pdf"), w=9,h=4)
ggsave(file=here("figs/jpeg","primary_composistion_box.jpeg"), w=9,h=4)







#############################################
#' ## Cell labels published in Elmentaite 2020 Dev. Cell
#############################################
## data imported in python (H5AD  file) and label exported to csv as seurat could not directly load the h5ad
dev_cell_labels<-read.csv(here("output","dev_cell_labels.csv"), sep="\t")

dev_cell_labels$index_edit<-gsub("-1","",dev_cell_labels$index)
dev_cell_labels_missing<-rownames(d10x.primary@meta.data)[which(!(rownames(d10x.primary@meta.data)%in%dev_cell_labels$index_edit))]
length(dev_cell_labels_missing)
head(dev_cell_labels_missing)

        ##simplify to one sample in both
        dev_cell_labels_T017<-dev_cell_labels[which(dev_cell_labels$Sample.name=="T017"),]
        dim(dev_cell_labels_T017)
        d10x.primary_T017<-subset(d10x.primary, subset = individual == "T017")
        dim(d10x.primary_T017@meta.data)
        
        ## all 1008 cells are in the raw data, 315 were filtered at strict MT threshold etc
        length(intersect(dev_cell_labels_T017$index_edit,rownames(d10x.primary_T017@meta.data)))

#samples in both
dev_cell_samples<-unique(dev_cell_labels$Sample.name)
primary_samples<-unique(d10x.primary@meta.data$individual)

intersect(dev_cell_samples,primary_samples)
length(intersect(dev_cell_samples,primary_samples))

dev_cell_samples[which(!(dev_cell_samples%in%primary_samples))]
primary_samples[which(!(primary_samples%in%dev_cell_samples))]

subset(d10x.primary, subset = individual %in% primary_samples[which(!(primary_samples%in%dev_cell_samples))])
# 17293 cells are from 7 samples not included in Elmentaite 2020 Dev. Cell
# therefore 26882/41771 cells were filtered in Elmentaite 2020 Dev. Cell but not so far in my data

## merge a dataframe for adding
dev_cell_labels_minimal<-dev_cell_labels[,c("index_edit","annotation_V2")]
dev_cell_labels_minimal_missing<-data.frame(index_edit=dev_cell_labels_missing, annotation_V2="unknown_filtered")
dev_cell_labels_minimal<-rbind(dev_cell_labels_minimal,dev_cell_labels_minimal_missing)

              
dev_cell_labels_minimal<-dev_cell_labels_minimal[match(colnames(d10x.primary), dev_cell_labels_minimal$index_edit),]
identical(colnames(d10x.primary), dev_cell_labels_minimal$index_edit)

d10x.primary<- AddMetaData(d10x.primary, dev_cell_labels_minimal$annotation_V2, col.name = "elmentaite_annotation_V2")


unique(d10x.primary@meta.data$elmentaite_annotation_V2)[order(unique(d10x.primary@meta.data$elmentaite_annotation_V2))]
which(unique(d10x.primary@meta.data$elmentaite_annotation_V2)[order(unique(d10x.primary@meta.data$elmentaite_annotation_V2))]=="unknown_filtered")

cols_manual<-c("#d7be3dff","#c8c380ff","#293163ff","#96bd7eff","#831d19ff","#4d6e7fff",
               "#19193bff","#476a78ff","#c8b189ff","#c8b189ff","#232a7aff","#e4ba8bff",
               "#c9c4a7ff","#946a20ff","#4c6d9aff","#4c6d9aff","#4c6d9aff","#821c17ff",
               "#6a5661ff","#e8c723ff","#4c6d9aff","#b43d1fff","#86451eff","#53729cff",
               "#68692eff","#c9c4a7ff","#c9c4a7ff","#941b1dff","#94171bff","#94692aff",
               "#242b7bff","#c8c380ff","#42665aff","#8e5486ff","#8a056bff","#c88ebcff",
               "#232a7aff","#75a59eff","#2d4456ff","#96abc6ff","lightgrey","#3a5826ff")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols = gg_color_hue(42)
cols[41] = "#d7be3dff"


DimPlot(d10x.primary, reduction = "umap", group.by = "elmentaite_annotation_V2", pt.size=0.25, label=T)+scale_color_manual(values=cols_manual)
ggsave(file=here("figs","primary_UMAP_cellLabel_devcell.pdf"), w=15,h=8)
ggsave(file=here("figs/jpeg","primary_UMAP_cellLabel_devcell.jpeg"), w=14,h=8)

DimPlot(d10x.primary, reduction = "umap", group.by = "elmentaite_annotation_V2", pt.size=0.25, label=F)+scale_color_manual(values=cols_manual)
ggsave(file=here("figs","primary_UMAP_cellLabel_devcell_nolab.pdf"), w=15,h=8)
ggsave(file=here("figs/jpeg","primary_UMAP_cellLabel_devcell_nolab.jpeg"), w=14,h=8)


DimPlot(d10x.primary, reduction = "umap", group.by = "compartment", pt.size=0.25)+
  scale_color_manual(values=c("#b2182b","#4393c3","#5aae61"))

DimPlot(d10x.primary, reduction = "umap", pt.size=0.25, label=T) ## saved above


## tsne
DimPlot(d10x.primary, reduction = "tsne", group.by = "elmentaite_annotation_V2", pt.size=0.25, label=T)+scale_color_manual(values=cols_manual)
DimPlot(d10x.primary, reduction = "tsne",  pt.size=0.25, label=T)

DimPlot(d10x.primary, reduction = "tsne", group.by = "compartment", pt.size=0.25)+
  scale_color_manual(values=c("#b2182b","#4393c3","#5aae61"))


#' highlight clusters
umap_mat<-as.data.frame(Embeddings(object = d10x.primary, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)
meta<-d10x.primary@meta.data
meta$cell<-rownames(meta)
plt_clust<-merge(meta, umap_mat,by="cell")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]}
cols = muted(gg_color_hue(31),l = 90, c = 30)
plt_clust<-plt_clust[!duplicated(plt_clust),]


# clust<-26
# clust<-12
# clust<-15
# 
# cols[clust+1] = gg_color_hue(31)[clust+1]

cols[26+1] = gg_color_hue(31)[26+1]
cols[12+1] = gg_color_hue(31)[12+1]
cols[15+1] = gg_color_hue(31)[15+1]


ggplot(plt_clust, aes(UMAP_1,UMAP_2, color=seurat_clusters))+
  geom_point(size=1.5, shape=20)+theme_bw()+
  theme_classic()+
  scale_color_manual(values=cols, name="Cluster")+th+theme(plot.title = element_text(size = 16, face = "bold"))





# ## how epi compartment lookin
# d10x.primary_epi<-subset(d10x.primary, subset = compartment == "epi")
# DimPlot(d10x.primary_epi, reduction = "umap", group.by = "elmentaite_annotation_V2", pt.size=0.25, label=T)+scale_color_manual(values=cols_manual)
# 
# #lynphatic endothelial seems confusing
# DotPlot(d10x.primary, features = epigenes<-c("EPCAM", "KRT8", "KRT18","COL1A1", "COL1A2", "COL6A1", "COL6A2", "VWF", "PLVAP", "CDH5", "S100B","CD52", "CD2", "CD3D", "CD3G", "CD3E", "CD79A", "CD79B", "CD14", "CD16", "CD68", "CD83", "CSF1R", "FCER1G"))


#' ## calling epithelial cell based n Developmental Cell calls
cell_label<-d10x.primary@meta.data[,c("seurat_clusters","elmentaite_annotation_V2")]
save(cell_label, file=here("output","cell_label_whole_clustering.RData"))



######### Add some rough compartments
d10x.primary@meta.data$general_type<-"Immune"
d10x.primary@meta.data$general_type[which(d10x.primary.raw@meta.data$seurat_wholePrimary_clusters%in%c(12,14,25))]<-"Epithelial"
d10x.primary@meta.data$general_type[which(d10x.primary.raw@meta.data$seurat_wholePrimary_clusters%in%c(5,10,13,15,16,22,23,28))]<-"Stromal"


DimPlot(d10x.primary, reduction = "pca", group.by = "general_type", pt.size=0.25, label=T)+scale_color_manual(values=c("#b2182b","#5aae61","#4393c3"))
ggsave(file=here("figs","primary_PCA_generalType.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","primary_PCA_generalType.jpeg"), w=6,h=4)


#################
#' # How do the POS and NEG look?
#################
d10x.primary@meta.data$POS_NEG<-"whole"
d10x.primary@meta.data$POS_NEG[which(d10x.primary@meta.data$individual%in%c("T110POS","T036POS"))]<-"POS"
d10x.primary@meta.data$POS_NEG[which(d10x.primary@meta.data$individual%in%c("T110NEG","T036NEG"))]<-"NEG"

DimPlot(d10x.primary, reduction = "umap", split = "POS_NEG", pt.size=0.25)
d10x.primary_POSNEG<-subset(d10x.primary, subset = POS_NEG %in% c("POS","NEG"))

DimPlot(d10x.primary_POSNEG, reduction = "umap", split = "individual",group.by = "general_type", pt.size=0.25)+ scale_color_manual(values=c("#b2182b","#5aae61","#4393c3"))
ggsave(file=here("figs","primary_pos_neg_UMAP.pdf"), w=20,h=6)
ggsave(file=here("figs/jpeg","primary_pos_neg_UMAP.jpeg"), w=20,h=6)

#seems like only 110 worked?
d10x.primary_T110<-subset(d10x.primary, subset = individual == "T110POS")

table(d10x.primary_T110@meta.data$elmentaite_annotation_V2)
table(d10x.primary@meta.data$general_type,d10x.primary@meta.data$individual)


# counts<-t(as.data.frame((tapply(d10x.primary@meta.data$orig.ident, list(d10x.primary@meta.data$general_type,d10x.primary@meta.data$individual), length))))
# counts[which(is.na(counts))]<-0
# counts<-(counts/rowSums(counts))*100
# counts<-as.data.frame(counts)
# counts$individual<-rownames(counts)
# counts<-melt(counts)
# 
# counts<-merge(counts, meta_primary, by.x="individual", by.y="donor")
# counts$Diagnosis[grep("POS",counts$individual)]<-"Positive"
# counts$Diagnosis[grep("NEG",counts$individual)]<-"Negative"
# 
# ggplot() + geom_bar(aes(y = value, x = individual, fill = variable), data = counts, color="black",stat="identity")+
#   facet_grid(~Diagnosis, scales="free_x",space="free_x")+theme_bw()+th+
#   scale_fill_manual(values=c("#b2182b","#5aae61","#4393c3"))+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# ggsave(file=here("figs","primary_composistion_dev_cell_refined.pdf"), w=10,h=4)
# ggsave(file=here("figs/jpeg","primary_composistion_dev_cell_refined.jpeg"), w=10,h=4)
# 
# ggplot(counts,aes(y = value, x = Diagnosis, fill = Diagnosis)) + geom_boxplot(outlier.shape=NA, color="black", fill="white")+
#   facet_wrap(~variable)+theme_bw()+th+scale_fill_manual(values=c("dodgerblue3","lightgrey","#5aae61","darkgrey","#b2182b","darkgoldenrod1"))+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+geom_point(shape=21, size=3)
# ggsave(file=here("figs","primary_composistion_dev_cell_refined_box.pdf"), w=7,h=4)
# ggsave(file=here("figs/jpeg","primary_composistion_dev_cell_refined_Box.jpeg"), w=7,h=4)
# 


#################
#### What cells are missing in Dev Cell paper and why those?
#################
# start with just the shared samples
d10x.primary_indevcell<-subset(d10x.primary, subset = individual %in% primary_samples[which((primary_samples%in%dev_cell_samples))])
