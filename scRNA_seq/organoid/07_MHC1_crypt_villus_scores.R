

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
d10x.stim<-readRDS(here("data","d10x_raw_merged.rds"))

## add cell type labels from split analysis
load(here("output","celltype_labels_organoid.RData"))
identical(rownames(cell_typelabel),rownames(d10x.stim@meta.data))
d10x.stim <- AddMetaData(d10x.stim, metadata = cell_typelabel)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.stim <- NormalizeData(d10x.stim,scale.factor = 10000, normalization.method = "LogNormalize")



d10x.stim <- AddModuleScore(
  object = d10x.stim,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score'
)

d10x.stim <- AddModuleScore(
  object = d10x.stim,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)

score_data<-d10x.stim@meta.data[,c("MHCI_score1","crypt_villis_score1","cluster_ID")]
head(score_data)

d10x.stim_norm<-readRDS(here("data","d10x_normalized.rds"))
d10x.stim_norm
identical(rownames(d10x.stim_norm@meta.data),rownames(d10x.stim@meta.data))
d10x.stim_norm <- AddMetaData(d10x.stim_norm, metadata = score_data)

#### data frame for plotting
umap_mat<-as.data.frame(Embeddings(object = d10x.stim_norm, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.stim_norm@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")
#cell_typelabel$cell<-rownames(cell_typelabel)
#plt<-merge(plt, cell_typelabel, by="cell")



#########
## MHCI
#########
FeaturePlot(d10x.stim_norm, features = "MHCI_score1",reduction = "umap", split.by ="orig.ident", min.cutoff = "q9", pt.size=1)

ggplot(plt, aes(UMAP_1,UMAP_2, color=MHCI_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))
ggsave(file=here("figs/jpeg","MHCI_score_UMAP.jpeg"), w=7, h=6)
ggsave(file=here("figs","MHCI_score_UMAP.pdf"), w=7, h=6)


ggplot(plt, aes(UMAP_1,UMAP_2, color=MHCI_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(cluster_ID~old.ident)
ggsave(file=here("figs/jpeg","MHCI_score_UMAP_facet.jpeg"), w=7, h=6)
ggsave(file=here("figs","MHCI_score_UMAP_facet.pdf"), w=7, h=6)


ggplot(plt, aes(cluster_ID,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=cluster_ID))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(individual~old.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))

ggplot(plt, aes(cluster_ID,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=cluster_ID))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~old.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))
# ggsave(file=here("figs/jpeg","MHCI_score_boxplot_celltype.jpeg"), w=12, h=5)
# ggsave(file=here("figs","MHCI_score_boxplot_celltype.pdf"), w=12, h=5)
# 
# 
ggplot(plt, aes(old.ident,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=old.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~cluster_ID)+scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))
# ggsave(file=here("figs/jpeg","MHCI_score_boxplot_stimulation.jpeg"), w=12, h=5)
# ggsave(file=here("figs","MHCI_score_boxplot_stimulation.pdf"), w=12, h=5)
# 


#######
### statistics
#######
library(stats)

# normality
model  <- lm(MHCI_score1 ~ orig.ident*cluster_ID, data = meta)
# qqplot(residuals(model))
# qqnorm(residuals(model), pch = 1, frame = FALSE)
# qqline(residuals(model), col = "steelblue", lwd = 2)

shapiro.test(residuals(model)[sample(1:length(residuals(model)),5000)])

## general treatment significance
print(summary(aov(MHCI_score1 ~ orig.ident, data = meta)))
pairwise<-pairwise.t.test(meta$MHCI_score1, meta$orig.ident,p.adjust.method = "BH", pool.sd = FALSE)
print(pairwise)
pairwise<-melt(pairwise$p.value)
pairwise[!is.na(pairwise$value),]


# Group the data by cell type do another anova
stats_plt_stim<-lapply(c("crypt","enterocyte","TA"), function(x){
  print("##################")
  cell_meta<-meta[which(meta$cluster_ID==x),]
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

comp_simple<-list(c("IFNg","NT"),c("NT","TNFa"),c("IFNg","TNFa"))

### CONFIRM P vale match the pairwise above
ggplot(plt, aes(old.ident,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=old.ident))+xlab("Treatment")+ylab("MHC I Score")+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~cluster_ID)+scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"), guide=F)+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#

ggsave(file=here("figs/jpeg","MHCI_score_boxplot_stimulation.jpeg"), w=12, h=5)
ggsave(file=here("figs","MHCI_score_boxplot_stimulation.pdf"), w=12, h=5)


# Group the data by stimulation do another anova
stats_plt<-lapply(c("TNFa","NT","IFNg"), function(x){
  print("##################")
  cell_meta<-meta[which(meta$orig.ident==x),]
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

comp_simple<-list(c("crypt","enterocyte"),c("crypt","TA"), c("enterocyte","TA"))


### CONFIRM P vale match the pairwise above
ggplot(plt, aes(cluster_ID,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=cluster_ID))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~old.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 2,  map_signif_level = T, color="grey60")#

ggsave(file=here("figs/jpeg","MHCI_score_boxplot_celltype.jpeg"), w=12, h=5)
ggsave(file=here("figs","MHCI_score_boxplot_celltype.pdf"), w=12, h=6)




# 
# #########
# ## Show between cell variability in score and IFN response
# #########
# range(d10x.stim@meta.data$MHCI_score1)
# d10x.stim_IFNg_lowMHC1_entero<-subset(d10x.stim, subset = orig.ident == "IFNg" & cluster_ID == "enterocyte" & MHCI_score1 <= -0.25)
# d10x.stim_IFNg_lowMHC1_entero
# d10x.stim_IFNg_highMHC1_entero<-subset(d10x.stim, subset = orig.ident == "IFNg" & cluster_ID == "enterocyte" & MHCI_score1 >= 1.6)
# d10x.stim_IFNg_highMHC1_entero
# 
# IFNg_lowMHC1_entero<-FetchData(object = d10x.stim_IFNg_lowMHC1_entero, vars = MHCI)
# IFNg_highMHC1_entero<-FetchData(object = d10x.stim_IFNg_highMHC1_entero, vars = MHCI)
# 
# IFNg_lowMHC1_entero$cell<-rownames(IFNg_lowMHC1_entero)
# IFNg_lowMHC1_entero<-melt(IFNg_lowMHC1_entero)
# IFNg_lowMHC1_entero$MHC1_score<-"low"
# 
# IFNg_highMHC1_entero$cell<-rownames(IFNg_highMHC1_entero)
# IFNg_highMHC1_entero<-melt(IFNg_highMHC1_entero)
# IFNg_highMHC1_entero$MHC1_score<-"high"
# 
# plt_sample_cells<-rbind(IFNg_lowMHC1_entero,IFNg_highMHC1_entero)
# plt_sample_cells$label<-paste("MHC I ",plt_sample_cells$MHC1_score, "\n(",plt_sample_cells$cell,")", sep="")
# 
# ggplot(plt_sample_cells, aes(variable, value))+geom_bar(stat="identity", fill="#238b45", color="black")+
#   theme_bw()+th+facet_wrap(~label)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# 
# 
# violin<-ggplot(plt[which(plt$cluster_ID=="enterocyte" & plt$old.ident=="IFNg"),], aes(1,MHCI_score1))+
#   geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,fill="#238b45")+
#   theme_bw()+th+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())+xlab("")+ylab("MHC I Score\nIFNg Enterocytes")
# 
# bar<-ggplot(plt_sample_cells, aes(variable, value))+geom_bar(stat="identity", fill="#238b45", color="black")+
#   theme_bw()+th+facet_wrap(~label, ncol=1)+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ylab("Count")
# 
# grid.arrange(violin, bar, ncol=2, widths=c(0.3,0.7))
# ggsave(file=here("figs/jpeg","MHCI_genes_high_low_enteroCells.jpeg"),grid.arrange(violin, bar, ncol=2, widths=c(0.3,0.7)), w=6, h=5)
# ggsave(file=here("figs","MHCI_genes_high_low_enteroCells.pdf"),grid.arrange(violin, bar, ncol=2, widths=c(0.3,0.7)), w=6, h=5)
# 





##################
## Crypt-Villus
##################
cor(plt$MHCI_score1, plt$crypt_villis_score1, method = "spearman")



# cell cycle and n feature regressed
FeaturePlot(d10x.stim_norm, features = "crypt_villis_score1",reduction = "umap", split.by ="orig.ident", min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.stim_norm, features = "crypt_villis_score1",reduction = "pca", split.by ="orig.ident", min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.stim_norm, features = "crypt_villis_score1",reduction = "pca",  min.cutoff = "q9", pt.size=1)

umap_mat<-as.data.frame(Embeddings(object = d10x.stim_norm, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.stim_norm@meta.data
meta$cell<-rownames(meta)
rm(d10x.stim_norm)
gc()

plt<-merge(meta, umap_mat, by="cell")
#cell_typelabel$cell<-rownames(cell_typelabel)
#plt<-merge(plt, cell_typelabel, by="cell")

ggplot(plt, aes(UMAP_1,UMAP_2, color=crypt_villis_score1))+
  geom_point(size=1.5)+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))
ggsave(file=here("figs/jpeg","crypt_villis_score_UMAP.jpeg"), w=7, h=6)
ggsave(file=here("figs","crypt_villis_score_UMAP.pdf"), w=7, h=6)



ggplot(plt, aes(cluster_ID,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=cluster_ID))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(individual~old.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))

ggplot(plt, aes(cluster_ID,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=cluster_ID))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~old.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))
# ggsave(file=here("figs/jpeg","crypt_villis_score_boxplot_celltype.jpeg"), w=12, h=5)
# ggsave(file=here("figs","crypt_villis_score_boxplot_celltype.pdf"), w=12, h=5)


ggplot(plt, aes(old.ident,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=old.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~cluster_ID)+scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))
# ggsave(file=here("figs/jpeg","crypt_villis_score_boxplot_stimulation.jpeg"), w=12, h=5)
# ggsave(file=here("figs","crypt_villis_score_boxplot_stimulation.pdf"), w=12, h=5)

# correlation MHC1 and cyrpt villus
ggplot(plt, aes(MHCI_score1, crypt_villis_score1))+geom_point(size=0.5)+facet_grid(orig.ident~cluster_ID)

gp<-plt %>% group_by(orig.ident,cluster_ID)
dplyr::summarize(gp, cor(MHCI_score1, crypt_villis_score1))


####
## Look at the normalized stimulation maps
####
IFNg.cells.integrated<-readRDS(here("data","IFNg_normalized_integrated.rds"))
IFNg.cells.integrated
IFNg.cells.integrated <- AddMetaData(IFNg.cells.integrated, metadata = score_data)


DimPlot(IFNg.cells.integrated, reduction = "umap", group.by = "cluster_ID", pt.size=0.5)+scale_color_manual(values=c("#6baed6","#238b45",
                                                                                                                     "#78c679","#f16913"))

FeaturePlot(IFNg.cells.integrated, features = "crypt_villis_score1", min.cutoff = "q9", pt.size=1)
FeaturePlot(IFNg.cells.integrated, features = "crypt_villis_score1", reduction="pca", min.cutoff = "q9", pt.size=1)



### pca raw

d10x.stim <- NormalizeData(d10x.stim)
d10x.stim <- FindVariableFeatures(d10x.stim, selection.method = "vst", nfeatures = 2000)
d10x.stim <- ScaleData(d10x.stim) #ScaleData(cells, vars.to.regress = c("nUMI","percent.mito","donor.id","S.Score","G2M.Score","batch_10X"))
d10x.stim <- RunPCA(d10x.stim, ndims.print = 1:10, nfeatures.print = 10)

FeaturePlot(d10x.stim, features = "crypt_villis_score1",reduction = "pca", min.cutoff = "q9", pt.size=1)
ggsave(file=here("figs","PCA_raw_crypt_villis_score1.pdf"), w=5.25,h=5)
ggsave(file=here("figs/jpeg","PCA_raw_crypt_villis_score1.jpeg"), w=5.25,h=5)




#######
### statistics
#######
library(stats)


# Group the data by cell type do another anova
stats_plt_stim<-lapply(c("crypt","enterocyte","TA"), function(x){
  print("##################")
  cell_meta<-meta[which(meta$cluster_ID==x),]
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

comp_simple<-list(c("IFNg","NT"),c("NT","TNFa"),c("IFNg","TNFa"))

### CONFIRM P vale match the pairwise above
ggplot(plt, aes(old.ident,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=old.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_grid(.~cluster_ID)+scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")#

ggsave(file=here("figs/jpeg","crypt_villis_score_boxplot_stimulation.jpeg"), w=12, h=5)
ggsave(file=here("figs","crypt_villis_score_boxplot_stimulation.pdf"), w=12, h=5)



# Group the data by stimulation do another anova
stats_plt<-lapply(c("TNFa","NT","IFNg"), function(x){
  print("##################")
  cell_meta<-meta[which(meta$orig.ident==x),]
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

comp_simple<-list(c("crypt","enterocyte"),c("crypt","TA"), c("enterocyte","TA"))


### CONFIRM P vale match the pairwise above
ggplot(plt, aes(cluster_ID,crypt_villis_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=cluster_ID))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~old.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#f16913"))+
  geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 2,  map_signif_level = T, color="grey60")#

ggsave(file=here("figs/jpeg","crypt_villis_score_boxplot_celltype.jpeg"), w=12, h=5)
ggsave(file=here("figs","crypt_villis_score_boxplot_celltype.pdf"), w=12, h=5)





#####################
## Just correlating no cell type grouping
#####################

# cell type colors
c("#6baed6","#238b45","#78c679","#f16913")
# treatment colors
c("cornflowerblue","grey80","firebrick4")



ggplot(plt, aes(crypt_villis_score1, MHCI_score1, fill=orig.ident, color=orig.ident))+
  geom_point(shape=21, color="black")+
  scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))+scale_color_manual(values=c("cornflowerblue","grey80","#cb2626"))+
  theme_bw()+th+
  stat_smooth(method="lm", se=F)

ggplot(plt, aes(crypt_villis_score1, MHCI_score1, fill=orig.ident))+geom_point(shape=21, color="white")+
  scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))+theme_bw()+th+
  stat_smooth(color="black")

ggplot(plt, aes(crypt_villis_score1, MHCI_score1, fill=cluster_ID))+geom_point(shape=21, color="white")+
  scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))+theme_bw()+th




scatter<-ggplot(plt, aes(crypt_villis_score1, MHCI_score1, fill=orig.ident, color=orig.ident))+geom_point(shape=21, color="black", size=1.5)+
  scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))+scale_color_manual(values=c("cornflowerblue","grey80","#cb2626"))+
  theme_bw()+th+theme(legend.position = "none",plot.margin = unit(c(0,0,1,2.45), "cm"))+
  stat_smooth(method="lm",se=F)


violin_celltype<-ggplot(plt, aes(reorder(cluster_ID, crypt_villis_score1), crypt_villis_score1, fill=cluster_ID))+
  geom_violin()+scale_x_discrete(position = "bottom") +coord_flip()+
  scale_fill_manual(values=c("#6baed6","#238b45","#f16913"))+theme_bw()+th+ylab("")+xlab("")+theme(axis.title.x=element_blank(),
                                                                                   axis.text.x=element_blank(),
                                                                                   axis.ticks.x=element_blank(),
                                                                                   legend.position = "none",
                                                                                   plot.margin = unit(c(1,0,0,1), "cm"))

box_diagnosis<-ggplot(plt, aes(orig.ident, MHCI_score1, fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))+theme_bw()+th+theme(axis.title.y =element_blank(),
                                          axis.text.y=element_blank(),
                                          axis.ticks.y=element_blank(),
                                          plot.margin = unit(c(0,0.1,1,0), "cm"))

mn_cryptvillus<-tapply(plt$crypt_villis_score1, plt$cluster_ID, mean)
cryptvillus_order<-names(mn_cryptvillus)[(order(mn_cryptvillus))]

plt_count<-melt(table(plt$orig.ident, plt$cluster_ID))
colnames(plt_count)<-c("orig.ident","cluster_ID","count")
plt_count$cluster_ID<-factor(plt_count$cluster_ID, cryptvillus_order)
placeholder<-ggplot(plt_count, aes(orig.ident, cluster_ID))+geom_tile(color="grey50", fill="white",size=0.4)+
  geom_text(aes(label=count),size=2.75,color="grey40")+ggtitle(" Cell Number")+
  theme_void()+theme(axis.title=element_blank(),
                     axis.text=element_blank(),
                     axis.ticks=element_blank(),
                     plot.margin = unit(c(0.4,2.8,0,0), "cm"),
                     plot.title = element_text(size = 10))
placeholder

grid_fig<-plot_grid(violin_celltype, placeholder, scatter, box_diagnosis, 
                    #align = 'hv',
                    ncol = 2, axis="lr", 
                    rel_heights = c(0.5,2),
                    rel_widths = c(4,1))

ggsave(file=here("figs", "MHCI_crpyt_summary_organoids.pdf"),grid_fig, w=10,h=10)
ggsave(file=here("figs/jpeg", "MHCI_crpyt_summary_organoids.jpeg"),grid_fig, w=10,h=10)



#### residuals
MHCI_CV_mod<-lm(plt$MHCI_score1~plt$crypt_villis_score1)

plt$resid<-residuals(summary(MHCI_CV_mod))


ggplot(plt, aes(crypt_villis_score1, resid, fill=orig.ident, color=orig.ident))+geom_point(shape=21, color="black", size=1.5)+
  scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))+theme_bw()+th

ggplot(plt, aes(orig.ident, resid, fill=orig.ident, color=orig.ident))+
  geom_boxplot(color="black")+
  scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))+theme_bw()+th

ggsave(file=here("figs", "MHCI_crpyt_residual_organoids.pdf"), w=5,h=5)
ggsave(file=here("figs/jpeg", "MHCI_crpyt_residual_organoids.jpeg"), w=5,h=5)



### boxplots split by cell type
plt$cluster_ID<-factor(plt$cluster_ID, (cryptvillus_order))

sig_df<-data.frame(cluster_ID=cryptvillus_order, sig=c("*","*","*"))
sig_df$cluster_ID<-factor(sig_df$cluster_ID, (cryptvillus_order))

ggplot(plt, aes(orig.ident,MHCI_score1))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.2,aes(fill=orig.ident))+
  theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  facet_wrap(~cluster_ID, nrow=1)+ scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))+
  geom_text(data=sig_df, aes(x=1.5, y=1.5, label=sig), size=8, color="grey40")
ggsave(file=here("figs", "MHCI_sortedby_crypt_villus_organoid.pdf"), w=6.5,h=3)
ggsave(file=here("figs/jpeg", "MHCI_sortedby_crypt_villus_organoid.jpeg"), w=6.5,h=3)







###############
#### IFNg receptors
###############
receptors<-FetchData(object = d10x.stim, vars = c("IFNGR1","IFNGR2"))
receptors$cell<-rownames(receptors)

plt<-merge(plt, receptors, by="cell")

scatter1<-ggplot(plt, aes(IFNGR1, MHCI_score1, fill=orig.ident, color=orig.ident))+geom_point(shape=21, color="black", size=1.5)+
  scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))+scale_color_manual(values=c("cornflowerblue","grey80","#cb2626"))+
  theme_bw()+th+stat_smooth(method="lm",se=F)+
  annotate("text", x=3, y=-0.8, label=paste("R =",round(cor(plt[which(plt$orig.ident=="IFNg"),"IFNGR1"],plt[which(plt$orig.ident=="IFNg"),"MHCI_score1"],method="spearman"),2)), color="cornflowerblue")+
  annotate("text", x=3, y=-0.9, label=paste("R =",round(cor(plt[which(plt$orig.ident=="NT"),"IFNGR1"],plt[which(plt$orig.ident=="NT"),"MHCI_score1"],method="spearman"),2)), color="grey80")+
  annotate("text", x=3, y=-1, label=paste("R =",round(cor(plt[which(plt$orig.ident=="TNFa"),"IFNGR1"],plt[which(plt$orig.ident=="TNFa"),"MHCI_score1"],method="spearman"),2)), color="#cb2626")

density1<-ggplot(plt, aes(IFNGR1))+geom_density(fill="#a6bddb")+
  theme_void()+th+xlab("IFNgR1 log(TP10K+1)")+theme(axis.title=element_blank(),
                                         axis.text=element_blank(),
                                         axis.ticks=element_blank(),
                                         legend.position = "none",
                                         panel.background = element_blank(),
                                         panel.grid.major = element_blank(), 
                                         panel.grid.minor = element_blank())

grid_fig1<-plot_grid(density1, scatter1,
                     align="v",
                    nrow = 2, axis="lr",
                    rel_heights = c(1,4))

grid_fig1


scatter2<-ggplot(plt, aes(IFNGR2, MHCI_score1, fill=orig.ident, color=orig.ident))+geom_point(shape=21, color="black", size=1.5)+
  scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))+scale_color_manual(values=c("cornflowerblue","grey80","#cb2626"))+
  theme_bw()+th+stat_smooth(method="lm",se=F)+
  annotate("text", x=3, y=-0.8, label=paste("R =",round(cor(plt[which(plt$orig.ident=="IFNg"),"IFNGR2"],plt[which(plt$orig.ident=="IFNg"),"MHCI_score1"],method="spearman"),2)), color="cornflowerblue")+
  annotate("text", x=3, y=-0.9, label=paste("R =",round(cor(plt[which(plt$orig.ident=="NT"),"IFNGR2"],plt[which(plt$orig.ident=="NT"),"MHCI_score1"],method="spearman"),2)), color="grey80")+
  annotate("text", x=3, y=-1, label=paste("R =",round(cor(plt[which(plt$orig.ident=="TNFa"),"IFNGR2"],plt[which(plt$orig.ident=="TNFa"),"MHCI_score1"],method="spearman"),2)), color="#cb2626")

density2<-ggplot(plt, aes(IFNGR1))+geom_density(fill="#a6bddb")+
  theme_void()+th+xlab("IFNGR2 log(TP10K+1)")+theme(axis.title=element_blank(),
                                                    axis.text=element_blank(),
                                                    axis.ticks=element_blank(),
                                                    legend.position = "none",
                                                    panel.background = element_blank(),
                                                    panel.grid.major = element_blank(), 
                                                    panel.grid.minor = element_blank())

grid_fig2<-plot_grid(density2, scatter2,
                     align="v",
                     nrow = 2, axis="lr",
                     rel_heights = c(1,4))

grid_fig2

grid.arrange(grid_fig1, grid_fig2, ncol=2)
ggsave(file=here("figs/jpeg","IFNg_receptors_organoids.jpeg"),grid.arrange(grid_fig1, grid_fig2, ncol=2), w=12, h=6)
ggsave(file=here("figs","IFNg_receptors_organoids.pdf"),grid.arrange(grid_fig1, grid_fig2, ncol=2), w=12, h=6)




## by cell type
ggplot(plt, aes(reorder(cluster_ID, crypt_villis_score1), IFNGR1, fill=cluster_ID))+
  geom_violin()+scale_x_discrete(position = "bottom") +
  scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))+
  theme_bw()+th+xlab("")

ggplot(plt, aes(reorder(cluster_ID, crypt_villis_score1), IFNGR2, fill=cluster_ID))+
  geom_violin()+scale_x_discrete(position = "bottom") +
  scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))+
  theme_bw()+th+xlab("")
