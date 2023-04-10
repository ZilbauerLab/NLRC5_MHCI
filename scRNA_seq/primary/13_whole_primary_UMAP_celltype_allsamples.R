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
library(RColorBrewer)
library(ggsignif)


options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))


d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))
print(d10x.primary)

# relabel meta data
d10x.primary@meta.data$orig.ident[which(d10x.primary@meta.data$orig.ident=="Ctrl")]<-"Control"
d10x.primary@meta.data$orig.ident[which(d10x.primary@meta.data$orig.ident=="neo")]<-"Neonatal"

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

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

d10x.primary <- CellCycleScoring(d10x.primary, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

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
print(d10x.primary)


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
cell_label$cluster_ID[which(cell_label$cluster_ID=="crypt")]<-"Stem"
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x.primary), cell_label$index),]
identical(colnames(d10x.primary), cell_label$index)


d10x.primary <- AddMetaData(d10x.primary, metadata = cell_label)
print(d10x.primary)

# ## rm unannotated
d10x.primary<-subset(d10x.primary, subset = cluster_ID != "unknown_filtered")

table(d10x.primary@meta.data$cluster_ID)


d10x.primary@meta.data$cell_type<-as.factor(d10x.primary@meta.data$cluster_ID)
d10x.primary@meta.data$cell_type<-factor(d10x.primary@meta.data$cell_type, levels = c(
  "BEST4 enterocyte","Stem","early enterocyte","enterocyte","enteroendocrine",
  "Goblet cell", "Paneth cell","Paneth (UC only)","Tuft",
  "Activated B cell","activated DC" ,"Activated T","Arterial endothelial cell" ,"B cell","CD4 T cell",
  "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
  "IgA plasma cell","IgG plasma cell","Lymphatic endothelial cell","Macrophage" ,"mast cells","Memory B cell" ,
  "Monocyte","pDC" ,"pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"TA","Tfh" ,
  "Treg", "Venous endothelial cell","Glial cell","myofibroblast","Mast cell"))

#' ### dev cell labels
cols_manual<-c(
  "#3D0AF2","#6521C1","#67D364","#367C34",
  "#B921C1","#308AC8",
  "#C86E30","#CFA218","#C12134","#238443","#810f7c" ,"#02818a","#8c510a" ,"#78c679","#67a9cf",
  "#3690c0","#8073ac","#b2abd2","#238443","#dd3497","#1d91c0","#addd8e","#014636",
  "#253494","#081d58","#bf812d","#7a0177" ,"#ce1256","#006d2c" ,
  "#a50f15","#542788" ,"#e31a1c","#e08214","#ef6548","#fdd49e" ,"#f46d43","#4575b4" ,
  "#016c59", "#bf812d","#c51b7d","#fb9a99","#fb6a4a")


names(cols_manual) <- levels(d10x.primary@meta.data$cell_type)
fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)




umap_mat<-as.data.frame(Embeddings(object = d10x.primary, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.primary@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")
print(table(plt$orig.ident))


####################################################################################################################

##############
## Expression all compartments
##############

MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")
crypt_villis = c("SEPP1", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8", "APOBEC1", "APOB", "APOA4", "APOA1", "NPC1L1", "EGFR", "KLF4", "ENPP3", "NT5E", "SLC28A2", "ADA")

## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))

## add cell type labels from split analysis
load(here("output","cell_label_whole_clustering.RData"))
identical(rownames(cell_label),rownames(d10x.primary@meta.data))
d10x.primary <- AddMetaData(d10x.primary, metadata = cell_label)

head(d10x.primary@meta.data)


##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.primary <- NormalizeData(d10x.primary,scale.factor = 10000, normalization.method = "LogNormalize")
d10x.primary_exp<-FetchData(object = d10x.primary, vars = c(MHCI,crypt_villis))


### save MHC I score too
d10x.primary <- AddModuleScore(
  object = d10x.primary,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score'
)

d10x.primary <- AddModuleScore(
  object = d10x.primary,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)

score_data<-d10x.primary@meta.data[,c("MHCI_score1","crypt_villis_score1")]



save(plt,d10x.primary_exp, score_data, file=here("data","Primary_UMAP_allsamples_plus_MHCI_raw.RData"))


############
## UMAP Score
############
load("/home/redgar/Documents/EBI/scRNAseq_codon/data/Primary_UMAP_allsamples_plus_MHCI_raw.RData")
plt$cell_type<-as.character(plt$cell_type)
plt$cell_type[which(plt$cluster_ID=="Paneth")]<-"Paneth"
plt[which(is.na(plt$cell_type)),]

d10x.primary_exp_MHC<-d10x.primary_exp[,which(colnames(d10x.primary_exp)%in%MHCI)]
d10x.primary_exp_MHC$cell<-rownames(d10x.primary_exp_MHC)
melt_exp<-melt(d10x.primary_exp_MHC)
plt_exp<-merge(plt, melt_exp, by="cell")

plt_exp_summary <- plt_exp %>% group_by(variable, cell_type, orig.ident) %>%
  dplyr::summarize(Mean = mean(value, na.rm=TRUE))
plt_exp_summary<- plt_exp_summary %>% group_by(variable) %>%  mutate(z_score = scale(Mean))
plt_exp_summary<-as.data.frame(plt_exp_summary)


epi<-c("BEST4 enterocyte","Stem","early enterocyte","enterocyte","enteroendocrine",
"Goblet cell", "Paneth cell","Paneth","Paneth (UC only)","Tuft","TA")
plt_exp_summary$compartment<-"Non-epithelial"
plt_exp_summary$compartment[which(plt_exp_summary$cell_type%in%epi)]<-"Epithelial"

plt_exp_summary$label<-paste(plt_exp_summary$orig.ident, plt_exp_summary$compartment, sep="\n")


## all cell types
ggplot(plt_exp_summary, aes(variable, cell_type, fill=z_score))+
  geom_tile()+facet_grid(label~., scales = "free_y", space = "free_y")+
  scale_fill_distiller(palette = "RdBu", name="Scaled\nMean\nExpression")+th+theme_classic()+
  ylab("")+xlab("")
ggsave(file="/home/redgar/Documents/EBI/scRNAseq_codon/figs/all_cells_MHCI.pdf", w=9, h=15)






##################
## select cell types
#################
## select cell types
plt_exp_select<-plt_exp[which(plt_exp$orig.ident %in% c("Control","CD")),]
plt_exp_select<-plt_exp_select[which(plt_exp_select$cell_type %in% c("activated DC","S1 fibroblasts", epi)),]
plt_exp_select<-plt_exp_select[which(!(plt_exp_select$cell_type == "Paneth (UC only)" )),]
plt_exp_select$cell_type<-factor(plt_exp_select$cell_type, levels=rev(c("Goblet cell", "enterocyte","BEST4 enterocyte","enteroendocrine","Tuft","early enterocyte",
                                                                        "TA","Paneth cell","Paneth","Stem","activated DC", "S1 fibroblasts")))

plt_exp_scaled_select <- plt_exp_select %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(value))
plt_exp_summary_select <- plt_exp_scaled_select %>% 
  group_by(variable, cell_type, orig.ident) %>%
  dplyr::summarize(Mean = mean(scaled, na.rm=TRUE))
plt_exp_summary_select<-as.data.frame(plt_exp_summary_select)

epi<-c("BEST4 enterocyte","Stem","early enterocyte","enterocyte","enteroendocrine",
       "Goblet cell", "Paneth cell","Paneth","Paneth (UC only)","Tuft","TA")
plt_exp_summary_select$compartment<-"Non-epithelial"
plt_exp_summary_select$compartment[which(plt_exp_summary_select$cell_type%in%epi)]<-"Epithelial"

plt_exp_summary_select$label<-paste(plt_exp_summary_select$orig.ident, plt_exp_summary_select$compartment, sep="\n")

cell_count<-plt_exp[,c("cell","orig.ident","cluster_ID")]
cell_count<-cell_count[!duplicated(cell_count),]
total_count<-tapply(cell_count$cell, list(cell_count$cluster_ID, cell_count$orig.ident), function(x) length(unique(x)))
total_count<-melt(total_count)

epi<-c("BEST4 enterocyte","Stem","early enterocyte","enterocyte","enteroendocrine",
       "Goblet cell", "Paneth cell","Paneth","Paneth (UC only)","Tuft","TA")
total_count$compartment<-"Non-epithelial"
total_count$compartment[which(total_count$Var1%in%epi)]<-"Epithelial"

total_count$label<-paste(total_count$Var2, total_count$compartment, sep="\n")

## select cell types
total_count_select<-total_count[which(total_count$Var2 %in% c("Control","CD")),]
total_count_select<-total_count_select[which(total_count_select$Var1 %in% c("activated DC","S1 fibroblasts", epi)),]
total_count_select<-total_count_select[which(!(total_count_select$Var1 == "Paneth (UC only)" )),]
total_count_select$cell_type<-factor(total_count_select$Var1, levels=rev(c("Goblet cell", "enterocyte","BEST4 enterocyte","enteroendocrine","Tuft","early enterocyte",
                                                                           "TA","Paneth cell","Paneth","Stem","activated DC", "S1 fibroblasts")))

## add mark for sig
load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/codon_july2022/redgar/scRNAseq_codon/data","primary_diff_genes.RData"))

diff_exp_sig<-diff_exp_all[which(diff_exp_all$p_val_adj<0.005),]

diff_exp_all_MHC<-diff_exp_sig[which(diff_exp_sig$gene%in%MHCI),]
diff_exp_all_MHC$celltype<-sapply(1:nrow(diff_exp_all_MHC), function(x) strsplit(diff_exp_all_MHC$cell.1[x],"_")[[1]][1])
diff_exp_all_MHC$celltype<-as.factor(diff_exp_all_MHC$celltype)
levels(diff_exp_all_MHC$celltype)[which(levels(diff_exp_all_MHC$celltype)%in%c("crypt"))]<-c("Stem")

plt_exp_summary_select$sig<-sapply(1:nrow(plt_exp_summary_select), function(x){
  cell_sig<-diff_exp_all_MHC[which(as.character(diff_exp_all_MHC$celltype)==as.character(plt_exp_summary_select$cell_type[x])),]
  if(nrow(cell_sig)==0){""}else{
    if(as.character(plt_exp_summary_select$variable[x])%in%cell_sig$gene){"*"}else{""}
  }
})


count_plt<-ggplot(total_count_select, aes(value, cell_type, fill=Var2, group=Var2))+
  geom_bar(stat = "identity", position=position_dodge(), color="black")+
  facet_grid(compartment~.,scales="free_y", space = "free_y")+th_present+theme_classic()+
  fillscale_diagnosis+xlab("Cell Count")+ylab("")+ theme(legend.position="top",
                                                         axis.title.y=element_blank(),
                                                         axis.text.y=element_blank(),
                                                         axis.ticks.y=element_blank(),
                                                         strip.background = element_blank(),
                                                         strip.text = element_blank(),
                                                         plot.margin = margin(1, 0.5, 0.25, 0, "cm"))

heat_plt<-ggplot(plt_exp_summary_select, aes(variable, cell_type, fill=Mean))+
  geom_tile(color="black")+facet_grid(compartment~orig.ident, scales = "free_y", space = "free_y")+
  th_present+theme_classic()+
  ylab("")+xlab("")+ theme(legend.position="top",
                           strip.background.y = element_blank(),
                           strip.text.y = element_blank())+
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev = T, name ="Mean\nScaled\nExpression",
                                   p3 = .1, p4 = .6, p1 = .2, p2 = .6) +
  geom_text(aes(label=sig))

plot_grid(heat_plt, count_plt, rel_widths = c(1,0.175), align = "h", axis = "tb")

ggsave(file="/home/redgar/Documents/EBI/scRNAseq_codon/figs/select_cells_MHCI_cellcount.pdf", w=16, h=5)



## non zero count
plt_exp_summary <- plt_exp %>% group_by(variable, cell_type, orig.ident) %>%
  dplyr::summarize(Mean = mean(value, na.rm=TRUE), cell_count = length(unique(cell)))
plt_exp_summary<- plt_exp_summary %>% group_by(variable) %>%  mutate(z_score = scale(Mean))
plt_exp_summary<-as.data.frame(plt_exp_summary)

epi<-c("BEST4 enterocyte","Stem","early enterocyte","enterocyte","enteroendocrine",
       "Goblet cell", "Paneth cell","Paneth","Paneth (UC only)","Tuft","TA")
plt_exp_summary$compartment<-"Non-epithelial"
plt_exp_summary$compartment[which(plt_exp_summary$cell_type%in%epi)]<-"Epithelial"

plt_exp_summary$label<-paste(plt_exp_summary$orig.ident, plt_exp_summary$compartment, sep="\n")


## select cell types
plt_exp_summary_select<-plt_exp_summary[which(plt_exp_summary$orig.ident %in% c("Control","CD")),]
plt_exp_summary_select<-plt_exp_summary_select[which(plt_exp_summary_select$cell_type %in% c("activated DC","S1 fibroblasts", epi)),]
plt_exp_summary_select<-plt_exp_summary_select[which(!(plt_exp_summary_select$cell_type == "Paneth (UC only)" )),]
plt_exp_summary_select$cell_type<-factor(plt_exp_summary_select$cell_type, levels=rev(c("Goblet cell", "enterocyte","BEST4 enterocyte","enteroendocrine","Tuft","early enterocyte",
                                                                                        "TA","Paneth cell","Paneth","Stem","activated DC", "S1 fibroblasts")))

plt_exp_summary_select<- plt_exp_summary_select %>% group_by(variable) %>%  mutate(z_score = scale(Mean))
plt_exp_summary_select<-as.data.frame(plt_exp_summary_select)

cell_count_nonzero<-plt_exp[which(plt_exp$value>0),]

plt_count_summary <- cell_count_nonzero %>% group_by(variable, cell_type, orig.ident) %>%
  dplyr::summarize(cell_count_gene=length(unique(cell)))
plt_count_summary<-as.data.frame(plt_count_summary)
plt_count_summary$compartment<-"Non-epithelial"
plt_count_summary$compartment[which(plt_count_summary$cell_type%in%epi)]<-"Epithelial"
plt_count_summary$label<-paste(plt_count_summary$orig.ident, plt_count_summary$compartment, sep="\n")


plt_exp_summary_select_count<-merge(plt_exp_summary_select,plt_count_summary[,c("cell_type","variable","orig.ident","cell_count_gene")], by=c("cell_type","variable","orig.ident"))
hist(plt_exp_summary_select_count$cell_count_gene)



ggplot(plt_exp_summary_select_count, aes(variable, cell_type, fill=z_score, size=cell_count_gene))+
  geom_point(shape=21)+facet_grid(label~., scales = "free_y", space = "free_y")+
  scale_fill_distiller(palette = "RdBu", name="Scaled\nMean\nExpression")+th+theme_classic()+
  ylab("")+xlab("")+scale_size_continuous(range = c(4,10), breaks = c(10,100,1000))

ggsave(file="/home/redgar/Documents/EBI/scRNAseq_codon/figs/select_cells_MHCI_dots_nonzero.pdf", w=14, h=6)




################
## NLRC5 Pseudo Bulk
################
epi<-c("BEST4 enterocyte","Stem","early enterocyte","enterocyte","enteroendocrine",
       "Goblet cell", "Paneth cell","Paneth","Paneth (UC only)","Tuft","TA")

plt_epi<-plt_exp[which(plt_exp$cell_type %in% epi),]

comp_simple<-list(c("CD","Control"))

ggplot(plt_epi, aes(orig.ident, value))+geom_boxplot()+geom_point()+
  facet_wrap(~variable, scales="free_y")+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#

ggplot(plt_epi, aes(orig.ident, log(value)))+geom_boxplot()+geom_point()+
  facet_wrap(~variable, scales="free_y")+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#


ggplot(plt_epi[which(plt_epi$variable=="NLRC5"),], aes(orig.ident, value))+geom_violin()+geom_boxplot()+geom_point()+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#

ggplot(plt_epi[which(plt_epi$variable=="NLRC5" & plt_epi$value>0),], aes(orig.ident, value))+geom_violin()+geom_boxplot(width=0.1)+geom_point()+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#


plt_epi_exp_summary <- plt_epi %>% group_by(variable, orig.ident) %>%
  dplyr::summarize(Mean = mean(value, na.rm=TRUE))
as.data.frame(plt_epi_exp_summary)


## number of zero
epi_NLRC5<-plt_epi[which(plt_epi$variable=="NLRC5"),]
table(epi_NLRC5$orig.ident)

epi_NLRC5_zero<-plt_epi[which(plt_epi$variable=="NLRC5" & plt_epi$value>0),]
table(epi_NLRC5_zero$orig.ident)

106/774
115/1465
18/357
