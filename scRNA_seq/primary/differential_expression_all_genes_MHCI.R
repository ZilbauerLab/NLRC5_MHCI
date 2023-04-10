
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
library(gtools)
library(ggsignif)


options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))



MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")


## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x.stim<-readRDS(here("data","d10x_raw_merged.rds"))

## add cell type labels from split analysis
load(here("output","celltype_labels_organoid.RData"))
identical(rownames(cell_typelabel),rownames(d10x.stim@meta.data))
d10x.stim <- AddMetaData(d10x.stim, metadata = cell_typelabel)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.stim <- NormalizeData(d10x.stim,scale.factor = 10000, normalization.method = "LogNormalize")



## testing factor
d10x.stim$cell_stim<-paste(d10x.stim$cluster_ID, d10x.stim$orig.ident, sep = "_")
Idents(d10x.stim) <- "cell_stim"

table(d10x.stim$cluster_ID, d10x.stim$orig.ident)

## just MHC I genes
d10x.stim<-subset(d10x.stim, features = MHCI)
d10x.stim

#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene, 
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). T

contrasts_celltype_stim<-rbind(
  combinations(n = 3, r = 2, v = unique(d10x.stim$cell_stim)[grep("crypt",unique(d10x.stim$cell_stim))], repeats.allowed = FALSE),
  combinations(n = 3, r = 2, v = unique(d10x.stim$cell_stim)[grep("TA",unique(d10x.stim$cell_stim))], repeats.allowed = FALSE),
  combinations(n = 3, r = 2, v = unique(d10x.stim$cell_stim)[grep("enterocyte",unique(d10x.stim$cell_stim))], repeats.allowed = FALSE))


# this is 193,941 tests across all comparisons (21549 genes, 9 comparisons)

diff_exp_all<-lapply(1:nrow(contrasts_celltype_stim), function(x){
  de<-FindMarkers(d10x.stim, ident.1 = contrasts_celltype_stim[x,1], ident.2 = contrasts_celltype_stim[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F, min.pct=0, logfc.threshold=0)
  print(paste(contrasts_celltype_stim[x,1],"vs", contrasts_celltype_stim[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_stim[x,1]
  de$cell.2<-contrasts_celltype_stim[x,2]
  de})


diff_exp_all_organoid<-do.call(rbind, diff_exp_all)






## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))

# no neonatal or posneg
d10x.primary<-subset(d10x.primary, subset = individual %in% c(
  "T017","T019","T176","T189","T197","T202","T203","T024","T036","T44","T057",
  "T160","T161","T175","T182","T184","T180"))
######
## add cell type labels from split analysis
######
#immune
load(here("output","immune_iterative_label.Rdata"))
#stromal
load(here("output","strom_celltype_label.Rdata"))
#epithelial
load(here("output","epi_celltype_label.Rdata"))

cell_label<-rbind(epi_cell_labels, immune_cell_labels, stromal_cell_labels)
cell_label$cluster_ID<-as.character(cell_label$cluster_ID)
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x.primary), cell_label$index),]
identical(colnames(d10x.primary), cell_label$index)

d10x.primary <- AddMetaData(d10x.primary, metadata = cell_label)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.primary <- NormalizeData(d10x.primary,scale.factor = 10000, normalization.method = "LogNormalize")

## remove doublets
d10x.primary<-subset(d10x.primary, subset = cluster_ID != "Doublet")
d10x.primary<-subset(d10x.primary, subset = cluster_ID != "Paneth (UC only)")
d10x.primary.cd.ctrl<-subset(d10x.primary, subset = orig.ident %in% c("CD","Ctrl"))

## just MHC I genes
d10x.primary.cd.ctrl<-subset(d10x.primary.cd.ctrl, features = MHCI)
d10x.primary.cd.ctrl

## testing factor
d10x.primary.cd.ctrl$cell_stim<-paste(d10x.primary.cd.ctrl$cluster_ID, d10x.primary.cd.ctrl$orig.ident, sep = "_")
Idents(d10x.primary.cd.ctrl) <- "cell_stim"

cell_types<-unique(d10x.primary.cd.ctrl$cluster_ID)

contrasts_celltype_stim<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x.primary.cd.ctrl$cell_stim[grep(cell_types[x],d10x.primary.cd.ctrl$cell_stim)], repeats.allowed = FALSE)}))

contrasts_celltype_stim[37,]<-c("enterocyte_CD","enterocyte_Ctrl")
contrasts_celltype_stim

diff_exp_all<-lapply(1:nrow(contrasts_celltype_stim), function(x){
  de<-FindMarkers(d10x.primary.cd.ctrl, ident.1 = contrasts_celltype_stim[x,1], ident.2 = contrasts_celltype_stim[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F, min.pct=0, logfc.threshold=0)
  if(nrow(de)>1){
  print(paste(contrasts_celltype_stim[x,1],"vs", contrasts_celltype_stim[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_stim[x,1]
  de$cell.2<-contrasts_celltype_stim[x,2]
  de}})


diff_exp_all_primary<-do.call(rbind, diff_exp_all)

save(diff_exp_all_primary, diff_exp_all_organoid, file=here("data","MHCI_only_diff_genes.RData"))


#load(file=here("../../../codon/scRNAseq_codon/data","MHCI_only_diff_genes.RData"))


lapply(MHCI,function(x) diff_exp_all_primary[grep(x,diff_exp_all_primary$gene),])
lapply(MHCI,function(x) diff_exp_all_organoid[grep(x,diff_exp_all_organoid$gene),])
lapply(MHCI,function(x) diff_exp_all[grep(x,diff_exp_all$gene),])

diff_exp_all[which(diff_exp_all$gene=="NLRC5"),]
diff_exp_all_organoid[which(diff_exp_all_organoid$gene=="NLRC5"),]
diff_exp_all_primary[which(diff_exp_all_primary$gene=="NLRC5"),]


diff_exp_all_primary[which(diff_exp_all_primary$cell.1=="enterocyte_CD"),]




##############
## Pseudo Bulk differential expression
##############