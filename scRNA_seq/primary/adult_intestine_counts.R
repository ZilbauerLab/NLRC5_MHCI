## Count table fro CIBERSORT Fliss bulk RNAseq

library(dplyr)
library(Seurat)
library(patchwork)
library(here)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(gtools)
library(colorspace)
library(ggsignif)
library(RColorBrewer)


options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))


adult_d10x.primary<-readRDS(here("data/adult_intestine","local.rds"))
adult_d10x.primary

adult_d10x.primary<-subset(adult_d10x.primary, subset = tissue %in% c("large intestine","rectum","small intestine"))
adult_d10x.primary<-subset(adult_d10x.primary, subset = Age_group %in% c("Adult"))
adult_d10x.primary

####################
## Count matrix
####################
# adult_d10x.primary <- adult_d10x.primary[, sample(colnames(adult_d10x.primary), size = 50000, replace=F)]
# adult_d10x.primary
# d10x_exp<-as.data.frame(adult_d10x.primary[["RNA"]]@data)

adult_d10x.primary <- FindVariableFeatures(adult_d10x.primary, nfeatures = 2000)
var_genes <- head(VariableFeatures(adult_d10x.primary), 2000)

d10x_exp <- GetAssayData(adult_d10x.primary)[var_genes,]

print(d10x_exp[1:5,1:5])
print(head(adult_d10x.primary))

adult_d10x.primary@meta.data$NAME<-rownames(adult_d10x.primary@meta.data)
adult_d10x.primary@meta.data<-adult_d10x.primary@meta.data[match(colnames(d10x_exp), adult_d10x.primary@meta.data$NAME),]

identical(adult_d10x.primary@meta.data$NAME, colnames(d10x_exp))

colnames(d10x_exp)<-adult_d10x.primary@meta.data$cell_type
d10x_exp_large<-d10x_exp[,which(adult_d10x.primary@meta.data$tissue=="large intestine")]
d10x_exp_small<-d10x_exp[,which(adult_d10x.primary@meta.data$tissue=="small intestine")]

save(d10x_exp_large, d10x_exp_small, file=here("data/adult_intestine","Adult_intestine_counts.RData"))

as.data.frame(as.matrix(d10x_exp_large))[1:5,1:5]

#load(here("../EBI/MHCI/more_single_cell/Adult_intestine_counts.RData"))
