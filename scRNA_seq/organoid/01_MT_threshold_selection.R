
#'---
#'title: scRNAseq seurat exploration
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---



#'### Load libraries
library(dplyr)
library(patchwork)
library(here)
library(ggplot2)
library(gridExtra)
library(cowplot)

library(Seurat)
library(reshape2)
#library(limma)


options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))

#################################
## Organoid
#################################

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

invisible(lapply(1:length(d10x.list), function(x){
  d10x.list[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list[[x]], pattern = "^MT-")}))

#' #### scrublet results
samples<-meta$SANGER.SAMPLE.ID

doublet_scores_all<-lapply(1:length(samples), function(x){
  ## scrublet output
  doublet_scores<-read.csv(paste(here("../ibd/output/"),"scrublet_output_table_",samples[x],".csv",sep=""))
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




lapply(c(10,20,30,40,50,90), function(mt){#10,20,30,40,50,90
    d10x.list.filtered<-d10x.list

    ## tinker with %MT level to see if there is a high MT cluster, once that dissappers that is the set point then the regression will handle the scattered MT high cells
    invisible(lapply(1:length(d10x.list.filtered), function(x){
      d10x.list.filtered[[x]] <<- subset(d10x.list.filtered[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < mt & predicted_doublet=="False")
    }))


    d10x.stim <- merge(d10x.list.filtered[[1]], d10x.list.filtered[[2]])
    d10x.stim <- merge(d10x.stim, d10x.list.filtered[[3]])

    d10x.stim <- NormalizeData(d10x.stim)
    d10x.stim <- FindVariableFeatures(d10x.stim, selection.method = "vst", nfeatures = 2000)
    d10x.stim <- ScaleData(d10x.stim) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))

    # dimension reduction
    d10x.stim <- RunPCA(d10x.stim, ndims.print = 1:10, nfeatures.print = 10)
    d10x.stim <- RunUMAP(d10x.stim, dims = 1:30)
    #d10x.stim <- RunTSNE(d10x.stim, dims = 1:30)

    ## Low quality cluster?
    ggsave(file=paste(here("figs/jpeg","organoid_MT"),mt,"_UMAP.jpeg", sep=""), grid.arrange(FeaturePlot(d10x.stim, features = "nCount_RNA", min.cutoff = "q9", pt.size=1),
                                                                                                   FeaturePlot(d10x.stim, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1),
                                                                                                   FeaturePlot(d10x.stim, features = "percent.mt", min.cutoff = "q9", pt.size=1)),w=6,h=15)

    })





####################################
## Primary
####################################
dataset_loc <- here("../ibd/data/raw/scRNAseq/scRNASEQ_DATA/Pediatric_samples/pediatric_cellranger302_GRCh38-3.0.0")
meta_primary<-read.csv(here("../ibd/data/raw/scRNAseq/scRNASEQ_DATA/Pediatric_samples/Pediatric_sample_tracker.csv"))
colnames(meta_primary)[2]<-"donor"

## From a TRIPP sheet Rasa shared
meta_primary_update<-read.csv(here("../ibd/data/primary_sc_TRIPP_update.csv"))
identical(meta_primary_update$Sample.name, meta_primary$donor)
meta_primary$Diagnosis<-c(meta_primary_update$Diagnosis_update)
meta_primary$Diagnosis[which(is.na(meta_primary$Diagnosis))]<-"neo"


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

invisible(lapply(1:length(d10x.list_primary), function(x){
  d10x.list_primary[[x]][["percent.mt"]] <<- PercentageFeatureSet(d10x.list_primary[[x]], pattern = "^MT-")}))


#' #### scrublet results
samples<-meta_primary$Sanger.Sample.ID

doublet_scores_all<-lapply(1:length(samples), function(x){
  ## scrublet output
  doublet_scores<-read.csv(paste(here("../ibd/output/"),"scrublet_output_table_",samples[x],".csv",sep=""))
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




lapply(c(10,20,30,40,50,90), function(mt){
  print(mt)
  d10x.list_primary_filtered<-d10x.list_primary
  
  invisible(lapply(1:length(d10x.list_primary_filtered), function(x){
    d10x.list_primary_filtered[[x]] <<- subset(d10x.list_primary_filtered[[x]], subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < mt & predicted_doublet=="False")
  }))
  
  d10x.list_primary_filtered_merge <- merge(d10x.list_primary_filtered[[1]], d10x.list_primary_filtered[[2]])
  d10x.list_primary_filtered_merge <- merge(d10x.list_primary_filtered_merge, d10x.list_primary_filtered[[3]])
  d10x.list_primary_filtered_merge <- merge(d10x.list_primary_filtered_merge, d10x.list_primary_filtered[[4]])
  
  d10x.list_primary_filtered_merge <- NormalizeData(d10x.list_primary_filtered_merge)
  d10x.list_primary_filtered_merge <- FindVariableFeatures(d10x.list_primary_filtered_merge, selection.method = "vst", nfeatures = 2000)
  d10x.list_primary_filtered_merge <- ScaleData(d10x.list_primary_filtered_merge) 
  
  d10x.list_primary_filtered_merge <- RunPCA(d10x.list_primary_filtered_merge, ndims.print = 1:10, nfeatures.print = 10)
  d10x.list_primary_filtered_merge <- RunUMAP(d10x.list_primary_filtered_merge, dims = 1:30)


  
  ggsave(file=paste(here("figs/jpeg","primary_MT"),mt,"_UMAP.jpeg", sep=""),  grid.arrange(FeaturePlot(d10x.list_primary_filtered_merge, features = "nCount_RNA", min.cutoff = "q9", pt.size=1),
                                                                                                   FeaturePlot(d10x.list_primary_filtered_merge, features = "nFeature_RNA", min.cutoff = "q9", pt.size=1),
                                                                                                   FeaturePlot(d10x.list_primary_filtered_merge, features = "percent.mt", min.cutoff = "q9", pt.size=1)) ,w=6,h=15)
})
  
