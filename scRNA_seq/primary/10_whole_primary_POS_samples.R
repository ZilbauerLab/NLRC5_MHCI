#'---
#'title: POS samples support for cell comp in DNAm IEC
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



d10x.primary<-readRDS(here("data","d10x_primary_normalized.rds"))


## add cell type labels from split analysis
#immune
load(here("output","immune_iterative_label.Rdata"))
#stromal
load(here("output","strom_celltype_label.Rdata"))
#epithelial
load(here("output","epi_celltype_label.Rdata"))

cell_label<-rbind(epi_cell_labels, immune_cell_labels, stromal_cell_labels)

cell_label$cluster_ID<-as.character(cell_label$cluster_ID)

cell_label$cluster_ID[which(cell_label$cluster_ID=="Neonatal_cell")]<-"Neonatal Epithelial"


cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x.primary), cell_label$index),]
identical(colnames(d10x.primary), cell_label$index)


d10x.primary <- AddMetaData(d10x.primary, metadata = cell_label)


# ## remove manually annotated doublets (180 cells)
d10x.primary<-subset(d10x.primary, subset = cluster_ID != "Doublet")

## add pos neg column
d10x.primary$pos_neg<-"whole"
d10x.primary$pos_neg[grep("POS",d10x.primary$individual)]<-"POS"
d10x.primary$pos_neg[grep("NEG",d10x.primary$individual)]<-"NEG"

## just POS NEG
d10x.primary_POSNEG<-subset(d10x.primary, subset = pos_neg != "whole")
table(d10x.primary_POSNEG$pos_neg,d10x.primary_POSNEG$individual)


#build a generic sampled control
d10x.primary_whole_controls<-subset(d10x.primary, subset = pos_neg == "whole" & orig.ident == "Control")
rm(d10x.primary)

## average proportions to report
count<-do.call(rbind, lapply(unique(d10x.primary_whole_controls@meta.data$individual), function(x){
  donor<-d10x.primary_whole_controls@meta.data[which(d10x.primary_whole_controls@meta.data$individual==x),]
  count<-as.data.frame(tapply(rownames(donor), list(donor$cluster_ID), length))
  count$donor<-x
  count$cell_type<-rownames(count)
  count
}))

colnames(count)[1]<-"count"

## average count of each cell type
types<-as.character(unique(count$cell_type))
types<-types[which(!(types%in%c("Paneth (UC only)","activated DC","pDC")))]# and average of on 1 makes a werid seurat object so excluding

cell_object_list<-lapply(1:length(types), function(x){
  print(types[x])
  n.cells<-round(mean(count$count[which(count$cell_type==types[x])], na.rm=T))
  print(n.cells)
  d10x.primary_whole_controls[, sample(colnames(d10x.primary_whole_controls)[which(d10x.primary_whole_controls$cluster_ID==types[x])], size = n.cells, replace=F)]})

length(cell_object_list)
print("cell object worked")

d10x.primary_whole_controls_average <- cell_object_list[[1]]
for (i in 2:length(x = cell_object_list)) {
  d10x.primary_whole_controls_average <- merge(d10x.primary_whole_controls_average, cell_object_list[[i]])
}
## 2422 cell in the control whole
d10x.primary_whole_controls_average
print("Average worked")


## get broad proportion to report
colnames(count)[1]<-"count"
count$cell_type<-factor(count$cell_type, levels = c("BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine", "Goblet cell", "Paneth cell","Tuft","TA",  "B cell","Activated B cell","Cycling B cell","FCER2 B cell","Memory B cell" ,
                                                    "Cycling plasma cell","IgA plasma cell","IgG plasma cell","Activated T","CD4 T cell", "CD8 T cell","gd T/NK cell", "Tfh" ,"Treg",  "Arterial endothelial cell", "Lymphatic endothelial cell", "pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"Venous endothelial cell","Glial cell","myofibroblast",
                                                    "Cycling myeloid cells","activated DC", "cDC1","cDC2", "Macrophage" ,"Monocyte", "pDC",  "Mast cell","unknown_filtered",   "Doublet","Neonatal CD4 T cell","Neonatal B cell","Neonatal Epithelial","Unidentified_only_UC"))
count$broad_cell_type<-count$cell_type
levels(count$broad_cell_type)<-c("Epi","Epi","Epi","Epi","Epi",   "Epi", "Epi","Epi","Epi", "B cell","B cell","B cell","B cell","B cell" , "plasma","plasma","plasma",  "T cell","T cell", "T cell","T cell", "T cell" ,"T cell",
                                 "Stromal", "Stromal", "Stromal","Stromal","Stromal","Stromal" ,"Stromal","Stromal","Stromal", "Myeloid","Myeloid", "Myeloid","Myeloid", "Myeloid" ,"Myeloid", "Myeloid", "Myeloid","unknown_filtered",
                                 "Doublet","T cell","B cell","Epi","Unidentified_only_UC")

count_broad_proportion<-do.call(rbind, lapply(unique(count$donor), function(x){
  donor<-count[which(count$donor==x),]
  count_broad<-as.data.frame(tapply(donor$count, list(donor$broad_cell_type), sum))
  count_broad$proportion<-count_broad$`tapply(donor$count, list(donor$broad_cell_type), sum)`/sum(count_broad$`tapply(donor$count, list(donor$broad_cell_type), sum)`, na.rm=T)
  count_broad$donor<-x
  count_broad$cell_type_broad<-rownames(count_broad)
  count_broad
}))

mean(count_broad_proportion$proportion[which(count_broad_proportion$cell_type_broad=="Epi")])
mean(count_broad_proportion$proportion[which(count_broad_proportion$cell_type_broad=="T cell")])
mean(count_broad_proportion$proportion[which(count_broad_proportion$cell_type_broad=="B cell")])
mean(count_broad_proportion$proportion[which(count_broad_proportion$cell_type_broad=="Stromal")])
mean(count_broad_proportion$proportion[which(count_broad_proportion$cell_type_broad=="plasma")])
mean(count_broad_proportion$proportion[which(count_broad_proportion$cell_type_broad=="Myeloid")])



#################
## JUST 110
#################
## just110
d10x.primary_POSNEG<-subset(d10x.primary_POSNEG, subset = individual %in% c("T110NEG","T110POS"))
table(d10x.primary_POSNEG$pos_neg,d10x.primary_POSNEG$individual)




## average sample combined with pos neg
table(d10x.primary_whole_controls_average$orig.ident, d10x.primary_whole_controls_average$individual)
d10x.primary_whole_controls_average$individual<-"average_cell_control"
table(d10x.primary_whole_controls_average$cluster_ID)



d10x.primary_POSNEG_onewhole <- merge(d10x.primary_POSNEG, d10x.primary_whole_controls_average)
VariableFeatures(d10x.primary_POSNEG_onewhole[["SCT"]]) <- rownames(d10x.primary_whole_controls[["SCT"]]@scale.data)
d10x.primary_POSNEG_onewhole <- ScaleData(d10x.primary_POSNEG_onewhole)
d10x.primary_POSNEG_onewhole <- RunPCA(d10x.primary_POSNEG_onewhole, ndims.print = 1:10, nfeatures.print = 10)
d10x.primary_POSNEG_onewhole <- RunUMAP(d10x.primary_POSNEG_onewhole, dims = 1:30)
d10x.primary_POSNEG_onewhole <- RunTSNE(d10x.primary_POSNEG_onewhole, dims = 1:30)



d10x.primary_POSNEG_onewhole@meta.data$cell_type<-as.factor(d10x.primary_POSNEG_onewhole@meta.data$cluster_ID)
d10x.primary_POSNEG_onewhole@meta.data$cell_type<-factor(d10x.primary_POSNEG_onewhole@meta.data$cell_type, levels = c("Unidentified_only_UC",
                                                                                                                      "BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine",
                                                                                                                      "Goblet cell", "Paneth cell","Tuft",
                                                                                                                      "Activated B cell","activated DC" ,"Activated T","Arterial endothelial cell" ,"B cell","CD4 T cell",
                                                                                                                      "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
                                                                                                                      "IgA plasma cell","IgG plasma cell","Lymphatic endothelial cell","Macrophage" ,"mast cells","Memory B cell" ,
                                                                                                                      "Monocyte","pDC" ,"pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"TA","Tfh" ,
                                                                                                                      "Treg", "Venous endothelial cell","Glial cell","myofibroblast","Mast cell","unknown_filtered","Doublet","Neonatal CD4 T cell","Neonatal B cell","Neonatal Epithelial"))

#' ### dev cell labels
cols_manual<-c(  "pink",
                 "#3D0AF2","#6521C1","#67D364","#367C34",
                 "#B921C1","#308AC8",
                 "#C86E30","#C12134","#238443","#810f7c" ,"#02818a","#8c510a" ,"#78c679","#67a9cf",
                 "#3690c0","#8073ac","#b2abd2","#238443","#dd3497","#1d91c0","#addd8e","#014636",
                 "#253494","#081d58","#bf812d","#7a0177" ,"#ce1256","#006d2c" ,
                 "#a50f15","#542788" ,"#e31a1c","#e08214","#ef6548","#fdd49e" ,"#f46d43","#4575b4" ,
                 "#016c59", "#bf812d","#c51b7d","#fb9a99","#fb6a4a","lightgrey","black","#c6dbef","#c7e9c0","pink")


names(cols_manual) <- levels(d10x.primary_POSNEG_onewhole@meta.data$cell_type)
fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)



# DimPlot(d10x.primary_POSNEG_onewhole, reduction = "umap", group.by = "cell_type",split.by = "individual", pt.size=0.25)+colcale_cols_manual
# 
# ggsave(file="../../figs/scRNAseq/primary/primary110_POSNEG_UMAP_celltype.pdf", w=16,h=4)
# ggsave(file="../../figs/scRNAseq/primary/jpeg/primary110_POSNEG_UMAP_celltype.jpeg", w=16,h=4)
# 
# 



count<-do.call(rbind, lapply(unique(d10x.primary_POSNEG_onewhole@meta.data$individual), function(x){
  donor<-d10x.primary_POSNEG_onewhole@meta.data[which(d10x.primary_POSNEG_onewhole@meta.data$individual==x),]
  count<-as.data.frame(tapply(rownames(donor), list(donor$cluster_ID), length))
  count$donor<-x
  count$cell_type<-rownames(count)
  count
}))

save(count, file=here("data","samp110_count_cell.RData"))

colnames(count)[1]<-"count"
count$cell_type<-factor(count$cell_type, levels = c("BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine",
                                                    "Goblet cell", "Paneth cell","Tuft","TA",
                                                    "B cell","Activated B cell","Cycling B cell","FCER2 B cell","Memory B cell" ,
                                                    "Cycling plasma cell","IgA plasma cell","IgG plasma cell",
                                                    "Activated T","CD4 T cell", "CD8 T cell","gd T/NK cell", "Tfh" ,"Treg",
                                                    "Arterial endothelial cell", "Lymphatic endothelial cell", "pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"Venous endothelial cell","Glial cell","myofibroblast",
                                                    "Cycling myeloid cells","activated DC", "cDC1","cDC2", "Macrophage" ,"Monocyte", "pDC",
                                                    "Mast cell"))


count$broad_cell_type<-count$cell_type
levels(count$broad_cell_type)<-c("Epi","Epi","Epi","Epi","Epi",
                                 "Epi", "Epi","Epi","Epi",
                                 "B cell","B cell","B cell","B cell","B cell" ,
                                 "plasma","plasma","plasma",
                                 "T cell","T cell", "T cell","T cell", "T cell" ,"T cell",
                                 "Stromal", "Stromal", "Stromal","Stromal","Stromal","Stromal" ,"Stromal","Stromal","Stromal",
                                 "Myeloid","Myeloid", "Myeloid","Myeloid", "Myeloid" ,"Myeloid", "Myeloid",
                                 "Myeloid")

ggplot(count, aes(cell_type, count, fill=cell_type))+geom_bar(stat = "identity", color="black")+facet_wrap(~donor)+fillscale_cols_manual+theme_bw()

table(count$broad_cell_type, count$cell_type)


ggplot(count, aes(cell_type, count, fill=cell_type))+
  geom_bar(stat = "identity", color="black")+facet_grid(broad_cell_type~donor)+fillscale_cols_manual+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")+xlab("")+ylab("Cell Count")


d10x.primary_POSNEG_onewhole$cell_type<-factor(d10x.primary_POSNEG_onewhole$cell_type, levels = c("BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine",
                                                                                                  "Goblet cell", "Paneth","Paneth (UC only)","Tuft","TA",
                                                                                                  "B cell","Activated B cell","Cycling B cell","FCER2 B cell","Memory B cell" ,
                                                                                                  "Cycling plasma cell","IgA plasma cell","IgG plasma cell",
                                                                                                  "Activated T","CD4 T cell", "CD8 T cell","gd T/NK cell", "Tfh" ,"Treg",
                                                                                                  "Arterial endothelial cell", "Lymphatic endothelial cell", "pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"Venous endothelial cell","Glial cell","myofibroblast",
                                                                                                  "Cycling myeloid cells","activated DC", "cDC1","cDC2", "Macrophage" ,"Monocyte", "pDC",
                                                                                                  "Mast cell"))

d10x.primary_POSNEG_onewhole$broad_cell_type<-d10x.primary_POSNEG_onewhole$cell_type


levels(d10x.primary_POSNEG_onewhole$broad_cell_type)<-c("Epithelial","Epithelial","Epithelial","Epithelial","Epithelial",
                                                        "Epithelial", "Epithelial","Epithelial","Epithelial","Epithelial",
                                                        "B cell","B cell","B cell","B cell","B cell" ,
                                                        "plasma","plasma","plasma",
                                                        "T cell","T cell", "T cell","T cell", "T cell" ,"T cell",
                                                        "Stromal", "Stromal", "Stromal","Stromal","Stromal","Stromal" ,"Stromal","Stromal","Stromal",
                                                        "Myeloid","Myeloid", "Myeloid","Myeloid", "Myeloid" ,"Myeloid", "Myeloid",
                                                        "Myeloid")

count_barplt<-d10x.primary_POSNEG_onewhole@meta.data

count_barplt$individual<-as.factor(count_barplt$individual)
table(count_barplt$individual)
levels(count_barplt$individual)<-c("Whole Biopsy (Control) n=2,454","EPCAM Negative n=2,416", "EPCAM Positive n=623")

save(count_barplt, file=here("data","barplot_count_cell.RData"))

#load(here("../../../codon/scRNAseq_codon/data","barplot_count_cell.RData"))
count_barplt$cell_type<-as.character(count_barplt$cell_type)
count_barplt[which(is.na(count_barplt$cell_type)),]$cell_type<-count_barplt[which(is.na(count_barplt$cell_type)),]$cluster_ID
count_barplt[which(count_barplt$cell_type=="crypt"),]$cell_type<-"Stem"
count_barplt[which(count_barplt$cell_type=="Paneth"),]$broad_cell_type<-"Epithelial"
count_barplt<-count_barplt[which(!(count_barplt$cell_type=="Paneth (UC only)")),]


#' ### dev cell labels

cols_manual<-c(  
  "#87d435ff","cornflowerblue","#238b43ff","#5e1e6fff",
  "#aa0044ff","#d9c026ff", "#CFA218","#f16916ff","#a13da1ff",
  "#238443","#810f7c" ,"#02818a","#8c510a" ,"#78c679","#67a9cf",
  "#3690c0","#8073ac","#b2abd2","#238443","#dd3497","#1d91c0","#addd8e","#014636",
  "#253494","#081d58","#bf812d","#7a0177" ,"#ce1256","#006d2c" ,
  "#a50f15","#542788" ,"#e31a1c","#e08214","#ef6548","#fdd49e" ,"#4575b4" ,
  "#016c59", "#bf812d","#c51b7d","#fb9a99","#fb6a4a")


names(cols_manual) <- c(
  "BEST4 enterocyte","Stem","enterocyte","enteroendocrine",
  "Goblet cell", "Paneth","Paneth (UC only)","TA","Tuft",
  "Activated B cell","activated DC" ,"Activated T","Arterial endothelial cell" ,"B cell","CD4 T cell",
  "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
  "IgA plasma cell","IgG plasma cell","Lymphatic endothelial cell","Macrophage" ,"mast cells","Memory B cell" ,
  "Monocyte","pDC" ,"pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"Tfh" ,
  "Treg", "Venous endothelial cell","Glial cell","myofibroblast","Mast cell")

fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
  
ggplot(count_barplt) +
  geom_bar(aes(x = broad_cell_type, fill = factor(cell_type)), color="black")+facet_wrap(~individual)+fillscale_cols_manual+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("Broad Cell Category")+ylab("Cell Count")+th

ggsave(file=here("figs","POSNEG110_proportion_celltype.pdf"), w=16,h=6)
ggsave(file=here("figs/jpeg","POSNEG110_proportion_celltype.jpeg"), w=16,h=6)



count_broad<-do.call(rbind, lapply(unique(count$donor), function(x){
  donor<-count[which(count$donor==x),]
  count_broad<-as.data.frame(tapply(donor$count, list(donor$broad_cell_type), sum))
  count_broad$donor<-x
  count_broad$cell_type_broad<-rownames(count_broad)
  count_broad
}))

count_broad

count_broad_proportion<-do.call(rbind, lapply(unique(count$donor), function(x){
  donor<-count[which(count$donor==x),]
  count_broad<-as.data.frame(tapply(donor$count, list(donor$broad_cell_type), sum))
  count_broad$proportion<-count_broad$`tapply(donor$count, list(donor$broad_cell_type), sum)`/sum(count_broad$`tapply(donor$count, list(donor$broad_cell_type), sum)`, na.rm=T)
  count_broad$donor<-x
  count_broad$cell_type_broad<-rownames(count_broad)
  count_broad
}))


#' 
#' 
#'               #################
#'               ## 36 and 110
#'               #################
#' 
#'               ## average sample combined with pos neg
#'               table(d10x.primary_whole_controls_average$orig.ident, d10x.primary_whole_controls_average$individual)
#'               d10x.primary_whole_controls_average$individual<-"average_cell_control"
#'               table(d10x.primary_whole_controls_average$cluster_ID)
#' 
#' 
#'               d10x.primary_POSNEG_onewhole <- merge(d10x.primary_POSNEG, d10x.primary_whole_controls_average)
#'               VariableFeatures(d10x.primary_POSNEG_onewhole[["SCT"]]) <- rownames(d10x.primary_POSNEG[["SCT"]]@scale.data)
#'               d10x.primary_POSNEG_onewhole <- ScaleData(d10x.primary_POSNEG_onewhole)
#'               d10x.primary_POSNEG_onewhole <- RunPCA(d10x.primary_POSNEG_onewhole, ndims.print = 1:10, nfeatures.print = 10)
#'               d10x.primary_POSNEG_onewhole <- RunUMAP(d10x.primary_POSNEG_onewhole, dims = 1:30)
#'               d10x.primary_POSNEG_onewhole <- RunTSNE(d10x.primary_POSNEG_onewhole, dims = 1:30)
#' 
#' 
#' 
#'               d10x.primary_POSNEG_onewhole@meta.data$cell_type<-as.factor(d10x.primary_POSNEG_onewhole@meta.data$cluster_ID)
#'               d10x.primary_POSNEG_onewhole@meta.data$cell_type<-factor(d10x.primary_POSNEG_onewhole@meta.data$cell_type, levels = c("Unidentified_only_UC",
#'                                                                                                                                     "BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine",
#'                                                                                                                                     "Goblet cell", "Paneth cell","Tuft",
#'                                                                                                                                     "Activated B cell","activated DC" ,"Activated T","Arterial endothelial cell" ,"B cell","CD4 T cell",
#'                                                                                                                                     "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
#'                                                                                                                                     "IgA plasma cell","IgG plasma cell","Lymphatic endothelial cell","Macrophage" ,"mast cells","Memory B cell" ,
#'                                                                                                                                     "Monocyte","pDC" ,"pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"TA","Tfh" ,
#'                                                                                                                                     "Treg", "Venous endothelial cell","Glial cell","myofibroblast","Mast cell","unknown_filtered","Doublet","Neonatal CD4 T cell","Neonatal B cell","Neonatal Epithelial"))
#' 
#'               #' ### dev cell labels
#'               cols_manual<-c(  "pink",
#'                                "#3D0AF2","#6521C1","#67D364","#367C34",
#'                                "#B921C1","#308AC8",
#'                                "#C86E30","#C12134","#238443","#810f7c" ,"#02818a","#8c510a" ,"#78c679","#67a9cf",
#'                                "#3690c0","#8073ac","#b2abd2","#238443","#dd3497","#1d91c0","#addd8e","#014636",
#'                                "#253494","#081d58","#bf812d","#7a0177" ,"#ce1256","#006d2c" ,
#'                                "#a50f15","#542788" ,"#e31a1c","#e08214","#ef6548","#fdd49e" ,"#f46d43","#4575b4" ,
#'                                "#016c59", "#bf812d","#c51b7d","#fb9a99","#fb6a4a","lightgrey","black","#c6dbef","#c7e9c0","pink")
#' 
#' 
#'               names(cols_manual) <- levels(d10x.primary_POSNEG_onewhole@meta.data$cell_type)
#'               fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
#'               colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)
#' 
#' 
#' 
#'               DimPlot(d10x.primary_POSNEG_onewhole, reduction = "umap", group.by = "cell_type",split.by = "individual", pt.size=0.25)+colcale_cols_manual
#' 
#'               ggsave(file="../../figs/scRNAseq/primary/primary_POSNEG_UMAP_celltype.pdf", w=20,h=4)
#'               ggsave(file="../../figs/scRNAseq/primary/jpeg/primary_POSNEG_UMAP_celltype.jpeg", w=20,h=4)
#' 
#' 
#' 
#' 
#'               count<-do.call(rbind, lapply(unique(d10x.primary_POSNEG_onewhole@meta.data$individual), function(x){
#'                 donor<-d10x.primary_POSNEG_onewhole@meta.data[which(d10x.primary_POSNEG_onewhole@meta.data$individual==x),]
#'                 count<-as.data.frame(tapply(rownames(donor), list(donor$cluster_ID), length))
#'                 count$donor<-x
#'                 count$cell_type<-rownames(count)
#'                 count
#'               }))
#' 
#' 
#' 
#'               colnames(count)[1]<-"count"
#'               count$cell_type<-factor(count$cell_type, levels = c("BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine",
#'                                                                   "Goblet cell", "Paneth cell","Tuft","TA",
#'                                                                   "B cell","Activated B cell","Cycling B cell","FCER2 B cell","Memory B cell" ,
#'                                                                   "Cycling plasma cell","IgA plasma cell","IgG plasma cell",
#'                                                                   "Activated T","CD4 T cell", "CD8 T cell","gd T/NK cell", "Tfh" ,"Treg",
#'                                                                   "Arterial endothelial cell", "Lymphatic endothelial cell", "pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"Venous endothelial cell","Glial cell","myofibroblast",
#'                                                                   "Cycling myeloid cells","activated DC", "cDC1","cDC2", "Macrophage" ,"Monocyte", "pDC",
#'                                                                   "Mast cell","unknown_filtered",
#'                                                                   "Doublet","Neonatal CD4 T cell","Neonatal B cell","Neonatal Epithelial","Unidentified_only_UC"))
#' 
#' 
#'               count$broad_cell_type<-count$cell_type
#'               levels(count$broad_cell_type)<-c("Epi","Epi","Epi","Epi","Epi",
#'                                                "Epi", "Epi","Epi","Epi",
#'                                                "B cell","B cell","B cell","B cell","B cell" ,
#'                                                "plasma","plasma","plasma",
#'                                                "T cell","T cell", "T cell","T cell", "T cell" ,"T cell",
#'                                                "Stromal", "Stromal", "Stromal","Stromal","Stromal","Stromal" ,"Stromal","Stromal","Stromal",
#'                                                "Myeloid","Myeloid", "Myeloid","Myeloid", "Myeloid" ,"Myeloid", "Myeloid",
#'                                                "Myeloid","unknown_filtered",
#'                                                "Doublet","T cell","B cell","Epithelial","Unidentified_only_UC")
#' 
#'               ggplot(count, aes(cell_type, count, fill=cell_type))+geom_bar(stat = "identity", color="black")+facet_wrap(~donor)+fillscale_cols_manual+theme_bw()
#' 
#'               table(count$broad_cell_type, count$cell_type)
#'               count<-count[which(count$broad_cell_type!="Unidentified_only_UC"),]
#'               count<-count[-grep("Neonatal", count$cell_type),]
#' 
#'               ggplot(count, aes(cell_type, count, fill=cell_type))+
#'                 geom_bar(stat = "identity", color="black")+facet_grid(broad_cell_type~donor)+fillscale_cols_manual+theme_bw()+
#'                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),legend.position = "none")+xlab("")+ylab("Cell Count")
#' 
#' 
#'               d10x.primary_POSNEG_onewhole$cell_type<-factor(d10x.primary_POSNEG_onewhole$cell_type, levels = c("BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine",
#'                                                                                                                 "Goblet cell", "Paneth cell","Tuft","TA",
#'                                                                                                                 "B cell","Activated B cell","Cycling B cell","FCER2 B cell","Memory B cell" ,
#'                                                                                                                 "Cycling plasma cell","IgA plasma cell","IgG plasma cell",
#'                                                                                                                 "Activated T","CD4 T cell", "CD8 T cell","gd T/NK cell", "Tfh" ,"Treg",
#'                                                                                                                 "Arterial endothelial cell", "Lymphatic endothelial cell", "pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"Venous endothelial cell","Glial cell","myofibroblast",
#'                                                                                                                 "Cycling myeloid cells","activated DC", "cDC1","cDC2", "Macrophage" ,"Monocyte", "pDC",
#'                                                                                                                 "Mast cell","unknown_filtered",
#'                                                                                                                 "Doublet","Neonatal CD4 T cell","Neonatal B cell","Neonatal Epithelial","Unidentified_only_UC"))
#' 
#'               d10x.primary_POSNEG_onewhole$broad_cell_type<-d10x.primary_POSNEG_onewhole$cell_type
#' 
#' 
#'               levels(d10x.primary_POSNEG_onewhole$broad_cell_type)<-c("Epi","Epi","Epi","Epi","Epi",
#'                                                                       "Epi", "Epi","Epi","Epi",
#'                                                                       "B cell","B cell","B cell","B cell","B cell" ,
#'                                                                       "plasma","plasma","plasma",
#'                                                                       "T cell","T cell", "T cell","T cell", "T cell" ,"T cell",
#'                                                                       "Stromal", "Stromal", "Stromal","Stromal","Stromal","Stromal" ,"Stromal","Stromal","Stromal",
#'                                                                       "Myeloid","Myeloid", "Myeloid","Myeloid", "Myeloid" ,"Myeloid", "Myeloid",
#'                                                                       "Myeloid","unknown_filtered",
#'                                                                       "Doublet","T cell","B cell","Epithelial","Unidentified_only_UC")
#' 
#'               count_barplt<-d10x.primary_POSNEG_onewhole@meta.data
#'               count_barplt<-count_barplt[which(count_barplt$broad_cell_type!="Unidentified_only_UC"),]
#'               count_barplt<-count_barplt[-grep("Neonatal", count_barplt$cell_type),]
#' 
#'               count_barplt$individual<-as.factor(count_barplt$individual)
#'               table(count_barplt$individual)
#'               levels(count_barplt$individual)<-c("Whole Biopsy (Control) n=2,347","EPCAM Negative T036 n=281", "EPCAM T036 Positive n=3,066", "EPCAM Negative T110 n=2,406", "EPCAM T110 Positive n=623")
#' 
#' 
#'               ggplot(count_barplt) +
#'                 geom_bar(aes(x = broad_cell_type, fill = factor(cell_type)), color="black")+facet_wrap(~individual)+fillscale_cols_manual+theme_bw()+
#'                 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+xlab("Broad Cell Category")+ylab("Cell Count")+th
#' 
#' 
#'               ggsave(file=here("figs","POSNEG_proportion_celltype.pdf"), w=16,h=10)
#'               ggsave(file=here("figs/jpeg","POSNEG_proportion_celltype.jpeg"), w=16,h=10)
#' 
#' 
#' 
#'               count_broad<-do.call(rbind, lapply(unique(count$donor), function(x){
#'                 donor<-count[which(count$donor==x),]
#'                 count_broad<-as.data.frame(tapply(donor$count, list(donor$broad_cell_type), sum))
#'                 count_broad$donor<-x
#'                 count_broad$cell_type_broad<-rownames(count_broad)
#'                 count_broad
#'               }))
#' 
#'               count_broad