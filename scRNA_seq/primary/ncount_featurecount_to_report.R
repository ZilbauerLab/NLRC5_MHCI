


#'### Load libraries
library(dplyr)
library(Seurat)
library(here)

options(stringsAsFactors = FALSE)

d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))

primary_QC_data<-d10x.primary@meta.data
primary_QC_data<-primary_QC_data[-grep("NEG|POS", primary_QC_data$individual),]
primary_QC_data<-primary_QC_data[which(primary_QC_data$orig.ident%in%c("CD","Ctrl","UC")),]

print("Primary nCount_RNA")
mean(primary_QC_data$nCount_RNA)
range(primary_QC_data$nCount_RNA)

print("Primary nFeature_RNA")
mean(primary_QC_data$nFeature_RNA)
range(primary_QC_data$nFeature_RNA)



d10x.stim<-readRDS(here("data","d10x_raw_merged.rds"))

organoid_QC_data<-d10x.primary@meta.data

print("Primary nCount_RNA")
mean(organoid_QC_data$nCount_RNA)
range(organoid_QC_data$nCount_RNA)

print("Primary nFeature_RNA")
mean(organoid_QC_data$nFeature_RNA)
range(organoid_QC_data$nFeature_RNA)
