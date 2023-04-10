

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
crypt_villis = c("SEPP1", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8", "APOBEC1", "APOB", "APOA4", "APOA1", "NPC1L1", "EGFR", "KLF4", "ENPP3", "NT5E", "SLC28A2", "ADA")


## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))

# no neonatal or posneg
d10x.primary<-subset(d10x.primary, subset = individual %in% c(
  "T017","T019","T176","T189","T197","T202","T203","T024","T036","T44","T057",
  "T160","T161","T175","T182","T184","T180"))


##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.primary <- NormalizeData(d10x.primary,scale.factor = 10000, normalization.method = "LogNormalize")



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


### load normalized for UMAP etc
d10x.stim_norm<-readRDS(here("data","d10x_primary_normalized.rds"))
# no neonatal or posneg
d10x.stim_norm<-subset(d10x.stim_norm, subset = individual %in% c(
  "T017","T019","T176","T189","T197","T202","T203","T024","T036","T44","T057",
  "T160","T161","T175","T182","T184","T180"))

identical(rownames(d10x.stim_norm@meta.data),rownames(d10x.primary@meta.data))
d10x.stim_norm <- AddMetaData(d10x.stim_norm, metadata = score_data)


## add cell type labels from split analysis
load(here("output","cell_label_whole_clustering.RData"))
cell_label<-cell_label[which(rownames(cell_label)%in%rownames(d10x.stim_norm@meta.data)),]
identical(rownames(cell_label),rownames(d10x.stim_norm@meta.data))
d10x.stim_norm <- AddMetaData(d10x.stim_norm, metadata = cell_label)

d10x.stim_norm@meta.data$general_type<-"Immune"
d10x.stim_norm@meta.data$general_type[which(d10x.stim_norm@meta.data$seurat_clusters%in%c(12,14,25))]<-"Epithelial"
d10x.stim_norm@meta.data$general_type[which(d10x.stim_norm@meta.data$seurat_clusters%in%c(5,10,13,15,16,22,23,28,24))]<-"Stromal"



DimPlot(d10x.stim_norm, reduction = "umap", group.by = "general_type", pt.size=0.25)+
  scale_color_manual(values=c("#b2182b","#4393c3","#5aae61"))






#########
## MHCI
#########
FeaturePlot(d10x.stim_norm, features = "MHCI_score1",reduction = "umap", min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.stim_norm, features = "MHCI_score1",reduction = "umap", split.by ="general_type", min.cutoff = "q9", pt.size=1)

umap_mat<-as.data.frame(Embeddings(object = d10x.stim_norm, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.stim_norm@meta.data
meta$cell<-rownames(meta)
rm(d10x.stim_norm)
gc()

plt<-merge(meta, umap_mat, by="cell")
#cell_typelabel$cell<-rownames(cell_typelabel)
#plt<-merge(plt, cell_typelabel, by="cell")

ggplot(plt, aes(UMAP_1,UMAP_2, color=MHCI_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))
ggsave(file=here("figs/jpeg","MHCI_score_UMAP_primary.jpeg"), w=7, h=6)
ggsave(file=here("figs","MHCI_score_UMAP_primary.pdf"), w=7, h=6)


ggplot(plt, aes(UMAP_1,UMAP_2, color=MHCI_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(general_type~old.ident)
ggsave(file=here("figs/jpeg","MHCI_score_UMAP_facet_primary.jpeg"), w=7, h=6)
ggsave(file=here("figs","MHCI_score_UMAP_facet_primary.pdf"), w=7, h=6)


ggplot(plt, aes(general_type,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=general_type))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~old.ident)+scale_fill_manual(values=c("#b2182b","#4393c3","#5aae61"))




# 
# 

ggplot(plt, aes(orig.ident,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=orig.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~general_type)+fillscale_diagnosis
# ggsave(file=here("figs/jpeg","MHCI_score_boxplot_stimulation_primary.jpeg"), w=12, h=5)
# ggsave(file=here("figs","MHCI_score_boxplot_stimulation_primary.pdf"), w=12, h=5)
# 


#######
### statistics
#######
library(stats)

# # normality
# model  <- lm(MHCI_score1 ~ orig.ident*general_type, data = meta)
# ggqqplot(residuals(model))
# qqnorm(residuals(model), pch = 1, frame = FALSE)
# qqline(residuals(model), col = "steelblue", lwd = 2)
# 
# shapiro.test(residuals(model)[sample(1:length(residuals(model)),5000)])



# Group the data by cell type do another anova
stats_plt_stim<-lapply(c("Immune","Epithelial","Stromal"), function(x){
  print("##################")
  cell_meta<-meta[which(meta$general_type==x),]
  print(unique(cell_meta$general_type))
  print(summary(aov(MHCI_score1 ~ orig.ident, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$MHCI_score1, cell_meta$orig.ident,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$orig.ident<-x
  pairwise
})

stats_plt_stim<-do.call(rbind, stats_plt_stim)
stats_plt_stim$Var1<-as.character(stats_plt_stim$Var1)
stats_plt_stim$Var2<-as.character(stats_plt_stim$Var2)

stats_plt_stim<-stats_plt_stim[-which(is.na(stats_plt_stim$value)),]

comp_simple<-list(c("CD","Control"),c("CD","UC"),c("Control","UC"))

### CONFIRM P vale match the pairwise above
ggplot(plt, aes(orig.ident,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=orig.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~general_type)+fillscale_diagnosis+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#

ggsave(file=here("figs/jpeg","MHCI_score_boxplot_stimulation_primary.jpeg"), w=12, h=5)
ggsave(file=here("figs","MHCI_score_boxplot_stimulation_primary.pdf"), w=12, h=5)

# Group the data by stimulation do another anova
stats_plt<-lapply(c("CD","Control","UC"), function(x){
  print("##################")
  cell_meta<-meta[which(meta$orig.ident==x),]
  print(unique(cell_meta$orig.ident))
  print(summary(aov(MHCI_score1 ~ general_type, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$MHCI_score1, cell_meta$general_type,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$orig.ident<-x
  pairwise
})

stats_plt<-do.call(rbind, stats_plt)
stats_plt$Var1<-as.character(stats_plt$Var1)
stats_plt$Var2<-as.character(stats_plt$Var2)

stats_plt<-stats_plt[-which(is.na(stats_plt$value)),]

comp_simple<-list(c("Epithelial","Immune"),c("Epithelial","Stromal"),c("Immune","Stromal"))


### CONFIRM P vale match the pairwise above
ggplot(plt, aes(general_type,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=general_type))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~old.ident)+scale_fill_manual(values=c("#b2182b","#4393c3","#5aae61"))+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.5,vjust = 0.5,
              textsize = 5,  map_signif_level = T, color="grey60")#

ggsave(file=here("figs/jpeg","MHCI_score_boxplot_celltype_primary.jpeg"), w=10, h=7)
ggsave(file=here("figs","MHCI_score_boxplot_celltype_primary.pdf"), w=10, h=7)


## simplier split just by compartment
pairwise.t.test(meta$MHCI_score1, meta$general_type,p.adjust.method = "BH", pool.sd = FALSE)

ggplot(plt, aes(general_type,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=general_type))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  scale_fill_manual(values=c("#b2182b","#4393c3","#5aae61"))+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#

ggsave(file=here("figs/jpeg","MHCI_score_boxplot_celltype_nosplit_primary.jpeg"), w=6, h=6)
ggsave(file=here("figs","MHCI_score_boxplot_celltype_nosplit_primary.pdf"), w=6, h=6)



###################
## Crypt villus
###################

ggplot(plt, aes(UMAP_1,UMAP_2, color=crypt_villis_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))
ggsave(file=here("figs/jpeg","crypt_villis_score_UMAP_primary.jpeg"), w=7, h=6)
ggsave(file=here("figs","crypt_villis_score_UMAP_primary.pdf"), w=7, h=6)


ggplot(plt, aes(UMAP_1,UMAP_2, color=crypt_villis_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(general_type~old.ident)
ggsave(file=here("figs/jpeg","crypt_villis_score_UMAP_facet_primary.jpeg"), w=7, h=6)
ggsave(file=here("figs","crypt_villis_score_UMAP_facet_primary.pdf"), w=7, h=6)


ggplot(plt, aes(general_type,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=general_type))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~old.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))
# ggsave(file="../../figs/scRNAseq/jpeg/crypt_villis_score_boxplot_celltype.jpeg"), w=12, h=5)
# ggsave(file="../../figs/scRNAseq/crypt_villis_score_boxplot_celltype.pdf"), w=12, h=5)
# 
# 
ggplot(plt, aes(orig.ident,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=orig.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~general_type)+fillscale_diagnosis
# ggsave(file=here("figs/jpeg","crypt_villis_score_boxplot_stimulation_primary.jpeg"), w=12, h=5)
# ggsave(file=here("figs","crypt_villis_score_boxplot_stimulation_primary.pdf"), w=12, h=5)
# 


#######
### statistics
#######
library(stats)

# # normality
# model  <- lm(crypt_villis_score1 ~ orig.ident*general_type, data = meta)
# ggqqplot(residuals(model))
# qqnorm(residuals(model), pch = 1, frame = FALSE)
# qqline(residuals(model), col = "steelblue", lwd = 2)
# 
# shapiro.test(residuals(model)[sample(1:length(residuals(model)),5000)])
# 


# Group the data by cell type do another anova
stats_plt_stim<-lapply(c("Immune","Epithelial","Stromal"), function(x){
  print("##################")
  cell_meta<-meta[which(meta$general_type==x),]
  print(unique(cell_meta$general_type))
  print(summary(aov(crypt_villis_score1 ~ orig.ident, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$crypt_villis_score1, cell_meta$orig.ident,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$orig.ident<-x
  pairwise
})

stats_plt_stim<-do.call(rbind, stats_plt_stim)
stats_plt_stim$Var1<-as.character(stats_plt_stim$Var1)
stats_plt_stim$Var2<-as.character(stats_plt_stim$Var2)

stats_plt_stim<-stats_plt_stim[-which(is.na(stats_plt_stim$value)),]

comp_simple<-list(c("CD","Control"),c("CD","UC"),c("Control","UC"))

### CONFIRM P vale match the pairwise above
ggplot(plt, aes(orig.ident,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=orig.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~general_type)+fillscale_diagnosis+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#

ggsave(file=here("figs/jpeg","crypt_villis_score_boxplot_stimulation_primary.jpeg"), w=12, h=5)
ggsave(file=here("figs","crypt_villis_score_boxplot_stimulation_primary.pdf"), w=12, h=5)

# Group the data by stimulation do another anova
stats_plt<-lapply(c("CD","Control","UC"), function(x){
  print("##################")
  cell_meta<-meta[which(meta$orig.ident==x),]
  print(unique(cell_meta$orig.ident))
  print(summary(aov(crypt_villis_score1 ~ general_type, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$crypt_villis_score1, cell_meta$general_type,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$orig.ident<-x
  pairwise
})

stats_plt<-do.call(rbind, stats_plt)
stats_plt$Var1<-as.character(stats_plt$Var1)
stats_plt$Var2<-as.character(stats_plt$Var2)

stats_plt<-stats_plt[-which(is.na(stats_plt$value)),]

comp_simple<-list(c("Epithelial","Immune"),c("Epithelial","Stromal"),c("Immune","Stromal"))


### CONFIRM P vale match the pairwise above
ggplot(plt, aes(general_type,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=general_type))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~old.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 2,  map_signif_level = T, color="grey60")#

ggsave(file=here("figs/jpeg","crypt_villis_score_boxplot_celltype_primary.jpeg"), w=12, h=5)
ggsave(file=here("figs","crypt_villis_score_boxplot_celltype_primary.pdf"), w=12, h=6)

















#################################
## Epithelial subtypes
#################################
## epithelial type labels
load(here("data","primary.epi.cells.RData"))

# no neonatal or posneg
primary.epi.cells<-subset(primary.epi.cells, subset = individual %in% c(
  "T017","T019","T176","T189","T197","T202","T203","T024","T036","T44","T057",
  "T160","T161","T175","T182","T184","T180"))

primary.epi.cells <- AddMetaData(primary.epi.cells, metadata = score_data)
#filter remaining immune (179 cells)
primary.epi.cells<-subset(primary.epi.cells, subset = cluster_ID != "B cell" & cluster_ID != "T cell")

# define colors
cols_manual_less<-c( "green","blue","#87d435ff",
                     "cornflowerblue",
                     "#238b43ff","#5e1e6fff","#aa0044ff",
                     "#d9c026ff","#CFA218",  
                     "#f16916ff","#a13da1ff")

DimPlot(primary.epi.cells, reduction="umap", group.by = "cluster_ID", pt.size=0.25)+scale_color_manual(values=cols_manual_less)
ggsave(file=here("figs/jpeg","UMAP_primary_epithelial_color_match.jpeg"), w=6, h=4)
ggsave(file=here("figs","UMAP_primary_epithelial_color_match.pdf"), w=6, h=4)

                                                                                                       
FeaturePlot(primary.epi.cells, features = "MHCI_score1",reduction = "umap", min.cutoff = "q9", pt.size=1)

umap_mat<-as.data.frame(Embeddings(object = primary.epi.cells, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-primary.epi.cells@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")

print("Correaltion to report")
cor(plt$MHCI_score1, plt$crypt_villis_score1, method = "spearman")


ggplot(plt, aes(UMAP_1,UMAP_2, color=MHCI_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))
ggsave(file=here("figs/jpeg","MHCI_score_UMAP_primary_epithelial.jpeg"), w=7, h=6)
ggsave(file=here("figs","MHCI_score_UMAP_primary_epithelial.pdf"), w=7, h=6)


ggplot(plt, aes(UMAP_1,UMAP_2, color=MHCI_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(cluster_ID~old.ident)
ggsave(file=here("figs/jpeg","MHCI_score_UMAP_facet_primary_epithelial.jpeg"), w=7, h=6)
ggsave(file=here("figs","MHCI_score_UMAP_facet_primary_epithelial.pdf"), w=7, h=6)




ggplot(plt, aes(cluster_ID,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.15,aes(fill=cluster_ID))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~old.ident)+scale_fill_manual(values=cols_manual_less)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#save this one cause too many test to plot
ggsave(file=here("figs/jpeg","MHCI_score_boxplot_celltype_primary_epithelial.jpeg"), w=12, h=8)
ggsave(file=here("figs","MHCI_score_boxplot_celltype_primary_epithelial.pdf"), w=12, h=8)
# 
# 
ggplot(plt, aes(orig.ident,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=orig.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~cluster_ID)+fillscale_diagnosis
# ggsave(file=here("figs/jpeg","MHCI_score_boxplot_stimulation_primary_epithelial.jpeg"), w=12, h=5)
# ggsave(file=here("figs","MHCI_score_boxplot_stimulation_primary_epithelial.pdf"), w=12, h=5)
# 


#######
### statistics
#######
library(stats)

# normality
model  <- lm(MHCI_score1 ~ orig.ident*cluster_ID, data = meta)
# ggqqplot(residuals(model))
# qqnorm(residuals(model), pch = 1, frame = FALSE)
# qqline(residuals(model), col = "steelblue", lwd = 2)

#shapiro.test(residuals(model)[sample(1:length(residuals(model)),5000)])

print("Main stat to report")
# Group the data by cell type do another anova
meta_noUConly<-meta[which(meta$cluster_ID!="Paneth (UC only)"),]

stats_plt_stim<-lapply(unique(meta_noUConly$cluster_ID), function(x){
  print("##################")
  cell_meta<-meta_noUConly[which(meta_noUConly$cluster_ID==x),]
  print(unique(cell_meta$cluster_ID))
  print(summary(aov(MHCI_score1 ~ orig.ident, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$MHCI_score1, cell_meta$orig.ident,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$orig.ident<-x
  pairwise
})

stats_plt_stim<-do.call(rbind, stats_plt_stim)
stats_plt_stim$Var1<-as.character(stats_plt_stim$Var1)
stats_plt_stim$Var2<-as.character(stats_plt_stim$Var2)

stats_plt_stim<-stats_plt_stim[-which(is.na(stats_plt_stim$value)),]
stats_plt_stim

comp_simple<-list(c("CD","Control"),c("CD","UC"),c("Control","UC"))

### CONFIRM P vale match the pairwise above
plt<-plt[which(!(plt$cluster_ID%in%c("Memory B cell","CD8 T cell","Paneth (UC only)"))),]
ggplot(plt, aes(orig.ident,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=orig.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~cluster_ID)+fillscale_diagnosis+
  geom_signif(comparisons = comp_simple, step_increase = 0.05,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#

ggsave(file=here("figs/jpeg","MHCI_score_boxplot_stimulation_primary_epithelial.jpeg"), w=12, h=10)
ggsave(file=here("figs","MHCI_score_boxplot_stimulation_primary_epithelial.pdf"), w=12, h=10)


## TOO MANY TEST TO PLOT


# Group the data by stimulation do another anova
stats_plt<-lapply(c("CD","Control","UC"), function(x){
  print("##################")
  if(x=="CD"){cell_meta<-meta_noUConly[which(meta_noUConly$orig.ident==x),]}else{
    if(x=="Control"){cell_meta<-meta_noUConly[which(meta_noUConly$orig.ident==x),]}else{
      if(x=="UC"){cell_meta<-meta_noUConly[which(meta_noUConly$orig.ident==x),]}else{
        if(x=="Neonatal"){cell_meta<-meta_noUConly[which(meta_noUConly$orig.ident==x & meta_noUConly$cluster_ID%in%c("crypt","Neonatal_cell")),]}}}}

  #cell_meta<-meta[which(meta$orig.ident==x),]
  print(unique(cell_meta$orig.ident))
  print(summary(aov(MHCI_score1 ~ cluster_ID, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$MHCI_score1, cell_meta$cluster_ID,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$orig.ident<-x
  pairwise
})

stats_plt<-do.call(rbind, stats_plt)
stats_plt$Var1<-as.character(stats_plt$Var1)
stats_plt$Var2<-as.character(stats_plt$Var2)

stats_plt<-stats_plt[-which(is.na(stats_plt$value)),]

stats_plt[which(stats_plt$value<0.005),]


## simplier split just by cell type
pairwise.t.test(meta_noUConly$MHCI_score1, meta_noUConly$cluster_ID,p.adjust.method = "BH", pool.sd = FALSE)

ggplot(plt, aes(cluster_ID,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=cluster_ID))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  scale_fill_manual(values=cols_manual_less)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file=here("figs/jpeg","MHCI_score_boxplot_celltype_nosplit_primary_epithelial.jpeg"), w=10, h=5)
ggsave(file=here("figs","MHCI_score_boxplot_celltype_nosplit_primary_epithelial.pdf"), w=10, h=5)

###################
## Crypt villus
###################


ggplot(plt, aes(UMAP_1,UMAP_2, color=crypt_villis_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))
ggsave(file=here("figs/jpeg","crypt_villis_score_UMAP_primary_epithelial.jpeg"), w=7, h=6)
ggsave(file=here("figs","crypt_villis_score_UMAP_primary_epithelial.pdf"), w=7, h=6)


ggplot(plt, aes(UMAP_1,UMAP_2, color=crypt_villis_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(cluster_ID~old.ident)
ggsave(file=here("figs/jpeg","crypt_villis_score_UMAP_facet_primary_epithelial.jpeg"), w=7, h=6)
ggsave(file=here("figs","crypt_villis_score_UMAP_facet_primary_epithelial.pdf"), w=7, h=6)


ggplot(plt, aes(cluster_ID,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.15,aes(fill=cluster_ID))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~old.ident)+scale_fill_manual(values=cols_manual_less)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
#save this one cause too many test to plot
ggsave(file=here("figs/jpeg","crypt_villis_score_boxplot_celltype_primary_epithelial.jpeg"), w=12, h=8)
ggsave(file=here("figs","crypt_villis_score_boxplot_celltype_primary_epithelial.pdf"), w=12, h=8)
# 
# 
ggplot(plt, aes(orig.ident,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=orig.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~cluster_ID)+fillscale_diagnosis
# ggsave(file=here("figs/jpeg","crypt_villis_score_boxplot_stimulation_primary_epithelial.jpeg"), w=12, h=5)
# ggsave(file=here("figs","crypt_villis_score_boxplot_stimulation_primary_epithelial.pdf"), w=12, h=5)
# 


#######
### statistics
#######

# Group the data by cell type do another anova
stats_plt_stim<-lapply(unique(meta_noUConly$cluster_ID), function(x){
  print("##################")
  cell_meta<-meta_noUConly[which(meta_noUConly$cluster_ID==x),]
  print(unique(cell_meta$cluster_ID))
  print(summary(aov(crypt_villis_score1 ~ orig.ident, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$crypt_villis_score1, cell_meta$orig.ident,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$orig.ident<-x
  pairwise
})

stats_plt_stim<-do.call(rbind, stats_plt_stim)
stats_plt_stim$Var1<-as.character(stats_plt_stim$Var1)
stats_plt_stim$Var2<-as.character(stats_plt_stim$Var2)

stats_plt_stim<-stats_plt_stim[-which(is.na(stats_plt_stim$value)),]

comp_simple<-list(c("CD","Control"),c("CD","UC"),c("Control","UC"))

### CONFIRM P vale match the pairwise above
ggplot(plt, aes(orig.ident,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=orig.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~cluster_ID)+fillscale_diagnosis+
  geom_signif(comparisons = comp_simple, step_increase = 0.05,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#

ggsave(file=here("figs/jpeg","crypt_villis_score_score_boxplot_stimulation_primary_epithelial.jpeg"), w=12, h=10)
ggsave(file=here("figs","crypt_villis_score_score_boxplot_stimulation_primary_epithelial.pdf"), w=12, h=10)


## TOO MANY TEST TO PLOT


# Group the data by stimulation do another anova
stats_plt<-lapply(c("CD","Control","UC"), function(x){
  print("##################")
  
  if(x=="CD"){cell_meta<-meta_noUConly[which(meta_noUConly$orig.ident==x),]}else{
    if(x=="Control"){cell_meta<-meta_noUConly[which(meta_noUConly$orig.ident==x),]}else{
      if(x=="UC"){cell_meta<-meta_noUConly[which(meta_noUConly$orig.ident==x),]}else{
        if(x=="Neonatal"){cell_meta<-meta_noUConly[which(meta_noUConly$orig.ident==x & meta_noUConly$cluster_ID%in%c("crypt","Neonatal_cell")),]}}}}
  
  #cell_meta<-meta[which(meta$orig.ident==x),]
  print(unique(cell_meta$orig.ident))
  print(summary(aov(crypt_villis_score1 ~ cluster_ID, data = cell_meta)))
  pairwise<-pairwise.t.test(cell_meta$crypt_villis_score1, cell_meta$cluster_ID,p.adjust.method = "BH", pool.sd = FALSE)
  print(pairwise)
  pairwise<-melt(pairwise$p.value)
  pairwise[!is.na(pairwise$value),]
  pairwise$orig.ident<-x
  pairwise
})

stats_plt<-do.call(rbind, stats_plt)
stats_plt$Var1<-as.character(stats_plt$Var1)
stats_plt$Var2<-as.character(stats_plt$Var2)

stats_plt<-stats_plt[-which(is.na(stats_plt$value)),]

stats_plt[which(stats_plt$value<0.005),]

## simplier split just by cell type
pairwise.t.test(meta_noUConly$crypt_villis_score1, meta_noUConly$cluster_ID,p.adjust.method = "BH", pool.sd = FALSE)

ggplot(plt, aes(cluster_ID,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=cluster_ID))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  scale_fill_manual(values=cols_manual_less)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(file=here("figs/jpeg","crypt_villis_score_score_boxplot_celltype_nosplit_primary_epithelial.jpeg"), w=10, h=5)
ggsave(file=here("figs","crypt_villis_score_score_boxplot_celltype_nosplit_primary_epithelial.pdf"), w=10, h=5)

