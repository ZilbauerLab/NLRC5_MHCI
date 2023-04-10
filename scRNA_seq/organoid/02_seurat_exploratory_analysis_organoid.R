
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
#' ORGANOIDS
##############################################
dataset_loc <- here("../ibd/data/raw/scRNAseq/scRNASEQ_DATA/Pediatric_organoids")

meta<-read.csv(here("../ibd/data/raw/scRNAseq/4918stdy_organoids.csv"))
meta$donor<-sapply(meta$SUPPLIER.SAMPLE.NAME, function(x) strsplit(x, " ")[[1]][1])
meta$passage<-sapply(meta$SUPPLIER.SAMPLE.NAME, function(x) strsplit(x, " ")[[1]][3])
meta$condition<-sapply(meta$SUPPLIER.SAMPLE.NAME, function(x) strsplit(x, " ")[[1]][4])


# Type on IFN vs INF
meta$SUPPLIER.SAMPLE.NAME<-gsub("INF","IFN", meta$SUPPLIER.SAMPLE.NAME)

#' Cohort table
meta[,c(17,18,19, 3)]

condition<-c("NT","IFNg","TNFa")

d10x.list<-lapply(1:length(condition),function(y) {
  ids <- as.character(meta$SANGER.SAMPLE.ID)[grep(condition[y],meta$SUPPLIER.SAMPLE.NAME)]
  
  d10x.data <- sapply(ids, function(i){
    d10x <- Read10X(file.path(dataset_loc,paste("out_",i, sep=""),"filtered_feature_bc_matrix"))
    colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
    d10x
  })
  
  d10x.data<-do.call("cbind",d10x.data)
  
  #' Initialize the Seurat object with the raw (non-normalized data).
  d10x <- CreateSeuratObject(counts = d10x.data, project = condition[y], min.cells = 3, min.features = 0)
  
  meta_cell<-data.frame(cell=colnames(d10x), individual=sapply(colnames(d10x), function(x) strsplit(x,"-")[[1]][2]))
  meta_cell_add<-merge(meta_cell, meta, by.x="individual", by.y="SANGER.SAMPLE.ID")
  #meta_cell_add<-meta_cell_add[match(meta_cell_add$cell, colnames(d10x)),]
  meta_cell_add<-meta_cell_add[match(colnames(d10x), meta_cell_add$cell),]
  
  identical(meta_cell_add$cell, colnames(d10x))
  
  d10x<- AddMetaData(d10x, meta_cell_add$donor, col.name = "individual")})

names(d10x.list)<-condition

## cell counts
plt_count_raw<-do.call(rbind,lapply(1:3, function(x) {
  df<-as.data.frame(tapply(rownames(d10x.list[[x]]@meta.data), list(d10x.list[[x]]@meta.data$individual, d10x.list[[x]]@meta.data$orig.ident), length))
  df$individual<-rownames(df)
  df$condition<-names(d10x.list)[x]
  colnames(df)<-c("cell_count","individual","condition")
  df}))
plt_count_raw

ggplot(plt_count_raw, aes(condition, cell_count))+
  geom_boxplot(fill="lightgrey")+geom_point()+
  theme_bw()+geom_text(aes(label=individual), hjust=-0.25)+ylim(0, 6000)+ylab("Total Cell Number")



#'## QC
#'The percentage of reads that map to the mitochondrial genome
#'Low-quality / dying cells often exhibit extensive mitochondrial contamination
#'We calculate mitochondrial QC metrics with the PercentageFeatureSet function, which calculates the percentage of counts originating from a set of features
#'We use the set of all genes starting with MT- as a set of mitochondrial genes

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list[[x]], pattern = "^MT-")}))

# Show QC metrics for the first 5 cells
head(d10x.list[[2]]@meta.data, 5)

#'Low-quality cells or empty droplets will often have very few genes
#'Cell doublets or multiplets may exhibit an aberrantly high gene count


# Visualize QC metrics
#nFeature number of unique genes
#nCount number of total molecules
plt_QC_data<-do.call(rbind, lapply(1:3, function(x) d10x.list[[x]]@meta.data))

ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Percent\nMitochondrial") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

ggsave(file=here("figs","ncount_nfeatures.pdf"), w=6,h=4)
ggsave(file=here("figs/jpeg","ncount_nfeatures.jpeg"), w=6,h=4)


ggplot(plt_QC_data,aes(percent.mt)) + geom_histogram(binwidth = 0.5) +
  geom_vline(xintercept = 10)+ theme_bw()+xlab("Percent Mitochondrial")+th

ggsave(file=here("figs","MT_features.pdf"), w=4,h=2)
ggsave(file=here("figs/jpeg","MT_features.jpeg"), w=4,h=2)


ggplot(plt_QC_data,aes(orig.ident,percent.mt)) + 
  geom_violin(fill="grey80", color="grey80") + geom_boxplot(aes(fill=orig.ident),width=0.1, outlier.size=0.25)+
  geom_hline(yintercept = 10)+ theme_bw()+ylab("Percent Mitochondrial")+th+xlab("Treatment")+
  scale_fill_manual(values=c("grey80","cornflowerblue","firebrick4"), guide=F)

ggsave(file=here("figs","MT_features_bytreatment.pdf"), w=4,h=3)
ggsave(file=here("figs/jpeg","MT_features_bytreatment.jpeg"), w=4,h=3)

summary(aov(plt_QC_data$percent.mt ~ plt_QC_data$orig.ident))



#' #### scrublet results
samples<-meta$SANGER.SAMPLE.ID

doublet_scores_all<-lapply(1:length(samples), function(x){
  ## scrublet output
  doublet_scores<-read.csv(paste(here("data/"),"scrublet_output_table_",samples[x],".csv",sep=""))
  doublet_scores$index_edit<-gsub("1",samples[x],doublet_scores$index)
  doublet_scores
})
doublet_scores_all<-do.call(rbind, doublet_scores_all)

## Match doublet information
doublet_scores_list<-lapply(1:length(d10x.list), function(x){
  doublet_scores_all<-doublet_scores_all[match(colnames(d10x.list[[x]]), doublet_scores_all$index_edit),]
  #identical(colnames(d10x.list[[x]]), doublet_scores_all$index_edit)
  doublet_scores_all
})

invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]]<<- AddMetaData(d10x.list[[x]], doublet_scores_list[[x]]$doublet_score, col.name = "doublet_score")
  d10x.list[[x]]<<- AddMetaData(d10x.list[[x]], doublet_scores_list[[x]]$predicted_doublet, col.name = "predicted_doublet")
}))

plt_QC_data<-do.call(rbind, lapply(1:3, function(x) d10x.list[[x]]@meta.data))

ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=doublet_score)) + 
  geom_point(size=0.75) + 
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Doublet\nScore") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

ggsave(file=here("figs","organoid_doublet_features.pdf"), w=5,h=3)
ggsave(file=here("figs/jpeg","organoid_doublet_features.jpeg"), w=5,h=3)

ggplot(plt_QC_data, aes(nCount_RNA,nFeature_RNA,colour=doublet_score)) + 
  geom_point(size=0.75) + facet_wrap(~predicted_doublet)+
  scale_color_gradientn(colors=c("black","blue","green2","red","yellow"),name="Doublet\nScore") +
  geom_hline(yintercept = 500) + xlab("Number of Total Molecules\n(nCount) ")+ylab("Number of Unique Genes\n(nFeature)")+
  geom_hline(yintercept = 6000) +theme_bw()+th

ggsave(file=here("figs","organoid_doublet_features_split.pdf"), w=12,h=4)
ggsave(file=here("figs/jpeg","organoid_doublet_features_spit.jpeg"), w=12,h=4)

table(plt_QC_data$predicted_doublet)
dim(plt_QC_data)
table(plt_QC_data$predicted_doublet, plt_QC_data$individual)







## what is 3 MAD from mean
median(plt_QC_data$nFeature_RNA)+(mad(plt_QC_data$nFeature_RNA)*3)
median(plt_QC_data$nFeature_RNA)-(mad(plt_QC_data$nFeature_RNA)*3)


#'We filter cells that have unique feature counts over 2,500 or less than 200
#'We filter cells that have >5% mitochondrial counts
d10x.list.raw<-d10x.list

invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]] <<- subset(d10x.list[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 & predicted_doublet=="False")
}))

d10x.list

## cell counts after QC
plt_count_QC<-do.call(rbind,lapply(1:3, function(x) {
  df<-as.data.frame(tapply(rownames(d10x.list[[x]]@meta.data), list(d10x.list[[x]]@meta.data$individual, d10x.list[[x]]@meta.data$orig.ident), length))
  df$individual<-rownames(df)
  df$condition<-names(d10x.list)[x]
  colnames(df)<-c("cell_count","individual","condition")
  df}))
plt_count_QC

#save to plot fancy
save(plt_count_raw, plt_count_QC, file=here("data","organoid_cell_count.RData"))

plt_count_raw$condition<-factor(plt_count_raw$condition, levels=c("NT","IFNg","TNFa"))
plt_count_QC$condition<-factor(plt_count_QC$condition, levels=c("NT","IFNg","TNFa"))

cell_count<-grid.arrange(ggplot(plt_count_raw, aes(condition, cell_count, fill=condition))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=individual), hjust=-0.25, size=3)+ylim(0, 6000)+xlab("Stimulation")+ylab("Total Cell Number")+th+
                           scale_fill_manual(values=c("grey80","cornflowerblue","firebrick4"), guide=F)+ggtitle("Before Quality Control"),
                         ggplot(plt_count_QC, aes(condition, cell_count, fill=condition))+
                           geom_boxplot()+geom_point()+
                           theme_bw()+geom_text(aes(label=individual), hjust=-0.25, size=3)+ylim(0, 6000)+xlab("Stimulation")+ylab("Total Cell Number")+th+
                           scale_fill_manual(values=c("grey80","cornflowerblue","firebrick4"), guide=F)+ggtitle("After Quality Control"), ncol=2)

ggsave(cell_count,file=here("figs","scCount.pdf"), w=8,h=4)
ggsave(cell_count,file=here("figs/jpeg","scCount.jpeg"), w=8,h=4)




# 
# $NT
# An object of class Seurat 
# 20617 features across 10992 samples within 1 assay 
# Active assay: RNA (20617 features, 0 variable features)
# 
# $IFNg
# An object of class Seurat 
# 19975 features across 5833 samples within 1 assay 
# Active assay: RNA (19975 features, 0 variable features)
# 
# $TNFa
# An object of class Seurat 
# 20299 features across 6362 samples within 1 assay 
# Active assay: RNA (20299 features, 0 variable features)

d10x.stim <- merge(d10x.list[[1]], d10x.list[[2]])
d10x.stim <- merge(d10x.stim, d10x.list[[3]])

# An object of class Seurat 
# 21549 features across 23187 samples within 1 assay 
# Active assay: RNA (21549 features, 0 variable features)

d10x.stim.raw<-d10x.stim

saveRDS(d10x.stim.raw, file = here("data","d10x_raw_merged.rds"))



d10x.stim<-readRDS(here("data","d10x_raw_merged.rds"))



d10x.stim <- NormalizeData(d10x.stim)
d10x.stim <- FindVariableFeatures(d10x.stim, selection.method = "vst", nfeatures = 2000)
d10x.stim <- ScaleData(d10x.stim) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

# dimension reduction
d10x.stim <- RunPCA(d10x.stim, ndims.print = 1:10, nfeatures.print = 10)
d10x.stim <- RunUMAP(d10x.stim, dims = 1:30)
d10x.stim <- RunTSNE(d10x.stim, dims = 1:30)



######################
## cell cycle gene expression
######################

## which PC associated to nfeature? 2
pca_mat<-as.data.frame(Embeddings(object = d10x.stim, reduction = "pca"))
pca_mat$cell<-rownames(pca_mat)
meta<-d10x.stim@meta.data
meta$cell<-rownames(meta)
plt<-merge(pca_mat, meta,by="cell")

sapply(2:51, function(x) cor(plt[,x],plt[,"nFeature_RNA"]))
FeaturePlot(d10x.stim, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
## PC2 tops genes include cell cycle genes: "TPX2"   "TOP2A"  "MKI67"  "CENPF"  "SMC4"   "TUBB4B"

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x.stim <- CellCycleScoring(d10x.stim, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

d10x.stim <- RunPCA(d10x.stim)
DimPlot(d10x.stim, reduction="pca")
ggsave(file=here("figs","cell_cycle_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","cell_cycle_PCA.jpeg"), w=5,h=5)

FeaturePlot(d10x.stim, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","nfeature_PCA.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","nfeature_PCA.jpeg"), w=5.25,h=5)


## regress out cell cycle and other covariates
#Transformed data will be available in the SCT assay, which is set as the default after running sctransform
#By default, sctransform accounts for cellular sequencing depth, or nUMIs.
d10x.stim <- suppressWarnings(SCTransform(d10x.stim, vars.to.regress = c("nFeature_RNA","S.Score", "G2M.Score"), verbose = FALSE))


# dimension reduction
d10x.stim <- RunPCA(d10x.stim, verbose = FALSE)
d10x.stim <- RunUMAP(d10x.stim, dims = 1:30)
d10x.stim <- RunTSNE(d10x.stim, dims = 1:30)


saveRDS(d10x.stim, file = here("data","d10x_normalized.rds"))


################################################################################################################


d10x.stim<-readRDS(here("data","d10x_normalized.rds"))


## visualize
DimPlot(d10x.stim, reduction = "umap", group.by = "Phase", pt.size=1)
DimPlot(d10x.stim, reduction = "umap", group.by = "Phase",split.by = "orig.ident", pt.size=1)


DimPlot(d10x.stim, reduction = "umap", group.by = "orig.ident", pt.size=0.5)+
  scale_color_manual(values=c("cornflowerblue","grey80","firebrick4"))
ggsave(file=here("figs","UMAP_raw_stimulation.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_raw_stimulation.jpeg"), w=7,h=6)

DimPlot(d10x.stim, reduction = "tsne", group.by = "orig.ident", pt.size=0.5)+
  scale_color_manual(values=c("cornflowerblue","grey80","firebrick4"))
ggsave(file=here("figs","tsne_raw_stimulation.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","tsne_raw_stimulation.jpeg"), w=7,h=6)

DimPlot(d10x.stim, reduction="pca", group.by = "orig.ident", pt.size=0.5)+
  scale_color_manual(values=c("cornflowerblue","grey80","firebrick4"))
ggsave(file=here("figs","organoid_PCA_stim.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","organoid_PCA_stim.jpeg"), w=5,h=4)

DimPlot(d10x.stim, reduction="pca", group.by = "orig.ident", pt.size=1)
DimPlot(d10x.stim, reduction = "umap", group.by = "individual", pt.size=1)

DimPlot(d10x.stim, reduction = "umap", split.by = "individual", group.by = "orig.ident", pt.size=1)+
  scale_color_manual(values=c("cornflowerblue","grey80","firebrick4"))
ggsave(file=here("figs","UMAP_raw_stimulation_individual.pdf"), w=5,h=4)
ggsave(file=here("figs/jpeg","UMAP_raw_stimulation_individual.jpeg"), w=21,h=6)

## Low quality cluster?
FeaturePlot(d10x.stim, features = "nCount_RNA", min.cutoff = "q9", pt.size=1)

FeaturePlot(d10x.stim, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","organoid_nfeature_umap_after_regression.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","organoid_nfeature_umap_after_regression.jpeg"), w=5.25,h=5)


FeaturePlot(d10x.stim, features = "doublet_score", pt.size=1)
ggsave(file=here("figs","organoid_doublet_score.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","organoid_doublet_score.jpeg"), w=5,h=5)

DimPlot(d10x.stim, reduction="pca", group.by="Phase")
ggsave(file=here("figs","organoid_cell_cycle_PCA_after_regression.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","organoid_cell_cycle_PCA_after_regression.jpeg"), w=5,h=5)

FeaturePlot(d10x.stim, features = "nFeature_RNA",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","organoid_nfeature_PCA_after_regression.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","organoid_nfeature_PCA_after_regression.jpeg"), w=5.25,h=5)



## genes in first pcs
pc_load<-Loadings(d10x.stim[['pca']])
print(d10x.stim[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(d10x.stim, dims = 1:2, reduction = "pca")



#########################
## Marker genes in Raw
#########################
genes<-unique(c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2",
                "PSMB9","PSMB8","B2M","MR1","CD1D","IRF1","NLRC5","TAP1","TAP2",
                "PSMB9","PSMB8","B2M","IRF1","NLRC5","SECTM1","CD320","CD3EAP",
                "CD177","CD74","CIITA","RELB",'HLA-DRA','HLA-DRB5','HLA-DRB1','HLA-DQA1','HLA-DQB1',
                'HLA-DQB1-AS1','HLA-DQA2', 'HLA-DQB2', 'HLA-DOB','HLA-DMB',
                'HLA-DMA','HLA-DOA','HLA-DPA1','HLA-DPB1',"LGR5","ASCL2","HELLS",
                "PCNA","DEFA6","MKI67","TOP2A","POU2F3","FCGBP","SOX4","FABP1","SLC2A2","SPIB","CHGA","SPINK4","KRT19"))

#log-normalized versions of corrected counts
d10x.stim.exp<-as.data.frame(d10x.stim[["SCT"]]@data)
d10x.stim.exp["FCGBP",1:10]

d10x.stim.exp.GOI<-d10x.stim.exp[genes,]
d10x.stim.exp.GOI$gene<-rownames(d10x.stim.exp.GOI)
d10x.stim.exp.GOI<-melt(d10x.stim.exp.GOI)#

umap_mat<-as.data.frame(Embeddings(object = d10x.stim, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.stim@meta.data
meta$cell<-rownames(meta)

plt<-merge(d10x.stim.exp.GOI, meta,by.x="variable", by.y="cell")
plt<-merge(plt, umap_mat,by.x="variable", by.y="cell")


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

grid.arrange(gene_plot("NLRC5"),gene_plot("TAP1"),gene_plot("PSMB9"),gene_plot("B2M"), ncol=2)

ggsave(grid.arrange(gene_plot("NLRC5"),gene_plot("TAP1"),gene_plot("PSMB9"),gene_plot("B2M"), ncol=2),
       file=here("figs","MHC1_raw_UMAP.pdf"), w=7,h=6)
ggsave(grid.arrange(gene_plot("NLRC5"),gene_plot("TAP1"),gene_plot("PSMB9"),gene_plot("B2M"), ncol=2),
       file=here("figs/jpeg","MHC1_raw_UMAP.jpeg"), w=7,h=6)











