

### Load libraries
library(dplyr)
library(Seurat)
library(here)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtools)
library(cowplot)
library(destiny)
library(viridis)


options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))

#' load full dataset
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))

## add cell type labels from split analysis
load(here("output","cell_label_whole_clustering.RData"))
identical(rownames(cell_label),rownames(d10x.primary@meta.data))
d10x.primary <- AddMetaData(d10x.primary, metadata = cell_label)

d10x.primary@meta.data$general_type<-"Immune"
d10x.primary@meta.data$general_type[which(d10x.primary@meta.data$seurat_clusters%in%c(12,14,25))]<-"Epithelial"
d10x.primary@meta.data$general_type[which(d10x.primary@meta.data$seurat_clusters%in%c(5,10,13,15,16,22,23,28,24))]<-"Stromal"

d10x.primary_epi_raw<-subset(d10x.primary, subset = general_type == "Epithelial")

load(here("data","primary.epi.cells.RData"))
cell_labels<-primary.epi.cells@meta.data

d10x.primary_epi_raw <- AddMetaData(d10x.primary_epi_raw, metadata = cell_labels)
#filter remaining immune (179 cells)
d10x.primary_epi_raw<-subset(d10x.primary_epi_raw, subset = cluster_ID != "B cell" & cluster_ID != "T cell"  & cluster_ID != "CD8 T cell" & cluster_ID != "Memory B cell")


##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.primary_epi_raw <- NormalizeData(d10x.primary_epi_raw,scale.factor = 10000, normalization.method = "LogNormalize")
d10x.primary_epi_raw

save(d10x.primary_epi_raw, file=here("data","d10x.primary_epi_raw.RData"))

load(here("data","d10x.primary_epi_raw.RData"))


########
## color pallette
########
cell_types<- c("Unidentified_only_UC",
               "BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine",
               "Goblet cell", "Paneth cell","Tuft",
               "Activated B cell","activated DC" ,"Activated T","Arterial endothelial cell" ,"B cell","CD4 T cell",
               "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
               "IgA plasma cell","IgG plasma cell","Lymphatic endothelial cell","Macrophage" ,"mast cells","Memory B cell" ,
               "Monocyte","pDC" ,"pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"TA","Tfh" ,
               "Treg", "Venous endothelial cell","Glial cell","myofibroblast","Mast cell","unknown_filtered","Doublet",
               "Neonatal CD4 T cell","Neonatal B cell","Neonatal Epithelial","Neonatal_cell")

cols_manual<-c(  "pink",
                 "#3D0AF2","#6521C1","#67D364","#367C34",
                 "#B921C1","#308AC8",
                 "#C86E30","#C12134","#238443","#810f7c" ,"#02818a","#8c510a" ,"#78c679","#67a9cf",
                 "#3690c0","#8073ac","#b2abd2","#238443","#dd3497","#1d91c0","#addd8e","#014636",
                 "#253494","#081d58","#bf812d","#7a0177" ,"#ce1256","#006d2c" ,
                 "#a50f15","#542788" ,"#e31a1c","#e08214","#ef6548","#fdd49e" ,"#f46d43","#4575b4" ,
                 "#016c59", "#bf812d","#c51b7d","#fb9a99","#fb6a4a","lightgrey","black","#c6dbef","#c7e9c0","pink","pink")


names(cols_manual) <- cell_types
cols_manual<-cols_manual[which(names(cols_manual)%in%unique(d10x.primary_epi_raw@meta.data$cluster_ID))]
fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)





############
## destiny
############
meta <- d10x.primary_epi_raw@meta.data
d10x.primary_epi_raw <- FindVariableFeatures(d10x.primary_epi_raw, selection.method = "vst", nfeatures = 3000)
d10x.primary_epi_raw.exp<-as.data.frame(d10x.primary_epi_raw[["RNA"]]@data)


## Main parameters
k=15
pcs=20
features_example=1000
topvargene <- head(VariableFeatures(d10x.primary_epi_raw), features_example)
d10x.primary_epi_raw.variable.exp<-d10x.primary_epi_raw.exp[which(rownames(d10x.primary_epi_raw.exp)%in%topvargene),]
dm<-DiffusionMap(t(d10x.primary_epi_raw.variable.exp), k=15, n_pcs=20)


############
#' Calculating actual DPT from the map
############
dpt<-DPT(dm)
plot(dpt)
ggsave(file=here("figs","DPT_primary_1000_15_20.pdf"),w=8, h=6)
ggsave(file=here("figs/jpeg","DPT_primary_1000_15_20.jpeg"),w=8, h=6)


plot(dpt, col_by="branch")
plot(dpt, root = 1, paths_to = c(2,3), col_by = 'branch', pch = 20)
plot(dpt, root = 1, paths_to = c(2,3), pch = 20)

plot(dpt, col_by = 'DPT1')
plot(dpt, col_by = 'DPT2')
plot(dpt, col_by = 'DPT3')


#' select the root cell as the tip of branch with stem cells
plt<-data.frame(DC1=eigenvectors(dm)[,1],
                DC2=eigenvectors(dm)[,2],
                DPT1=dpt$DPT1)

ggplot(plt, aes(DC1, DC2, color=DPT1))+geom_point()+theme_bw()

tips <- random_root(dm)
pt_vec <- dpt[tips[[1]], ]

dpt@tips[which(dpt@tips[,1]==T),]
dpt@branch[which(rownames(dpt@branch)%in%rownames(dpt@tips[which(dpt@tips[,1]==T),])),]

which(rownames(dpt@branch)%in%rownames(dpt@tips[which(dpt@tips[,1]==T),]))
dpt[453,453]

meta[which(rownames(dpt@branch)%in%rownames(dpt@tips[which(dpt@tips[,1]==T),])),]

######
## center of stem cluster
######
stem_center<-data.frame(DC1=eigenvectors(dm)[,1], DC2=eigenvectors(dm)[,2])
stem_center$cellname<-rownames(stem_center)
meta$cellname<-rownames(meta)
stem_center_celltype<-merge(stem_center, meta, by="cellname")
stem_center_celltype<-stem_center_celltype[which(stem_center_celltype$cluster_ID=="crypt"),]

mean(stem_center_celltype$DC1)
mean(stem_center_celltype$DC2)

stem_center_celltype[which(stem_center_celltype$DC1>(mean(stem_center_celltype$DC1)-0.0001) & 
                             stem_center_celltype$DC1<(mean(stem_center_celltype$DC1)+0.0001) & 
                             stem_center_celltype$DC2>(mean(stem_center_celltype$DC2)-0.0001) & 
                             stem_center_celltype$DC2<(mean(stem_center_celltype$DC2)+0.0001) ),]
#CATGACACATATGAGA-4918STDY7447825
which(rownames(dpt@branch)=="CATGACACATATGAGA-4918STDY7447825")

plot(dpt, col_by = 'DPT1147')
plot(dpt, root = 1)

plt<-data.frame(DC1=eigenvectors(dm)[,1],
                DC2=eigenvectors(dm)[,2],
                DPT=dpt$DPT1147)

plt$cellname<-rownames(plt)
meta$cellname<-rownames(meta)
plt_celltype<-merge(plt, meta, by="cellname")
save(plt_celltype, file=here("output","DPT_primary_result.RData"))

###########
## DPT crypt-villus
###########
load(here("output","DPT_primary_result.RData"))

#####
##score data
#####
MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")
crypt_villis = c("SEPP1", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8", "APOBEC1", "APOB", "APOA4", "APOA1", "NPC1L1", "EGFR", "KLF4", "ENPP3", "NT5E", "SLC28A2", "ADA")
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))
##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.primary <- NormalizeData(d10x.primary,scale.factor = 10000, normalization.method = "LogNormalize")

d10x.primary <- AddModuleScore(
  object = d10x.primary,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score')

d10x.primary <- AddModuleScore(
  object = d10x.primary,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score')

score_data<-d10x.primary@meta.data[,c("MHCI_score1","crypt_villis_score1")]

score_data$cellname<-rownames(score_data)
plt_celltype<-merge(plt_celltype,score_data,by="cellname")

#### no neo
plt_celltype$cluster_ID<-as.character(plt_celltype$cluster_ID)
plt_celltype$cluster_ID<-as.factor(plt_celltype$cluster_ID)
levels(plt_celltype$cluster_ID)<-c("BEST4\nEnterocyte","Stem","Enterocyte","Enteroendocrine","Goblet","Paneth","Paneth (UC only)","TA","Tuft")



########
## color pallette
########
cell_types<- c("BEST4\nEnterocyte","Stem","Early Enterocyte","Enterocyte","Enteroendocrine","Goblet",
               "Paneth","Paneth (UC only)","TA","Tuft")

cols_manual<-c( "#87d435ff",
                "cornflowerblue","#6ec58bff",
                "#238b43ff","#5e1e6fff","#aa0044ff",
                "#d9c026ff","#CFA218",  
                "#f16916ff","#a13da1ff")


names(cols_manual) <- cell_types
cols_manual<-cols_manual[which(names(cols_manual)%in%unique(plt_celltype$cluster_ID))]
fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)

## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}




cell_color<-ggplot()+ geom_point(aes(DC1,DC2, color=cluster_ID), plt_celltype)+
  theme_bw()+th_present+colcale_cols_manual+xlab("Diffusion Component 1")+ylab("Diffusion Component 2")
leg_umap_cell = get_leg(cell_color)
cell_color<- grid.arrange(cell_color+theme(legend.position = "none"),
                          leg_umap_cell, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("figs","DPT_primary_cellcolored_1000_15_20.pdf"),cell_color,w=8, h=6)
ggsave(file=here("figs/jpeg","DPT_primary_cellcolored_1000_15_20.jpeg"),cell_color,w=8, h=6)

DPT_color<-ggplot()+ geom_point(aes(DC1,DC2, color=DPT), plt_celltype)+
  theme_bw()+th_present+scale_color_viridis(option = "D")+xlab("Diffusion Component 1")+ylab("Diffusion Component 2")
leg_umap_DPT = get_leg(DPT_color)
DPT_color<- grid.arrange(DPT_color+theme(legend.position = "none"),
                         leg_umap_DPT, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("figs","DPT_primary_DPTcolored_1000_15_20.pdf"),DPT_color,w=8, h=6)
ggsave(file=here("figs/jpeg","DPT_primary_DPTcolored_1000_15_20.jpeg"),DPT_color,w=8, h=6)

ggplot()+ geom_point(aes(DC1,DC2, color=crypt_villis_score1), plt_celltype)+
  theme_bw()
ggplot()+ geom_point(aes(DC1,DC2, color=MHCI_score1), plt_celltype)+
  theme_bw()

plt_celltype$pseudotime_diffusionmap <- rank(plt_celltype$DPT)

ggplot(plt_celltype, aes(reorder(cluster_ID,pseudotime_diffusionmap, median ),pseudotime_diffusionmap, colour = cluster_ID)) +
  geom_point()+theme_bw()+  theme(axis.text.x = element_text(angle = 90, hjust=1))+colcale_cols_manual
ggsave(file=here("figs","DPT_primary_celltype.pdf"),DPT_color,w=8, h=6)
ggsave(file=here("figs/jpeg","DPT_primary_celltype.jpeg"),DPT_color,w=8, h=6)

ggplot(plt_celltype, aes(reorder(cluster_ID,pseudotime_diffusionmap, median ),crypt_villis_score1, colour = cluster_ID)) +
  geom_point()+theme_bw()+  theme(axis.text.x = element_text(angle = 90, hjust=1))+colcale_cols_manual
ggplot(plt_celltype, aes(reorder(cluster_ID,pseudotime_diffusionmap, median ),MHCI_score1, colour = cluster_ID)) +
  geom_point()+theme_bw()+  theme(axis.text.x = element_text(angle = 90, hjust=1))+colcale_cols_manual

ggplot(plt_celltype, aes(pseudotime_diffusionmap, crypt_villis_score1,colour = cluster_ID)) +
  geom_point()+theme_bw()+th_present+colcale_cols_manual+stat_smooth(method='loess', color="black")+
  xlab("Pseudotime")+ylab("Crypt-Villus Score")
ggsave(file=here("figs","DPT_primary_cryptcillus.pdf"),w=8, h=6)
ggsave(file=here("figs/jpeg","DPT_primary_cryptcillus.jpeg"),w=8, h=6)
ggplot(plt_celltype, aes(pseudotime_diffusionmap, MHCI_score1,colour = cluster_ID)) +
  geom_point()+theme_bw()+th_present+colcale_cols_manual+stat_smooth(method='loess', color="black")+
  xlab("Pseudotime")+ylab("MHC I Score")
ggsave(file=here("figs","DPT_primary_MHC1.pdf"),w=8, h=6)
ggsave(file=here("figs/jpeg","DPT_primary_MHC1.jpeg"),w=8, h=6)


### color by marker gene expression
markers<-FetchData(object = d10x.primary, vars = c("DEF6A","LGR5","FABP1","PLA2G2A"))
markers$cell<-rownames(markers)

plt_celltype<-merge(plt_celltype, markers,by.x="cellname", by.y="cell")

ggplot(plt_celltype, aes(pseudotime_diffusionmap, LGR5)) +geom_point(aes(colour = cluster_ID))+theme_bw()+th+colcale_cols_manual+stat_smooth(method='loess', color="black")
ggsave(file=here("figs","DPT_primary_LGR5.pdf"),w=8, h=6)
ggsave(file=here("figs/jpeg","DPT_primary_LGR5.jpeg"),w=8, h=6)
ggplot(plt_celltype, aes(pseudotime_diffusionmap, PLA2G2A)) +geom_point(aes(colour = cluster_ID))+theme_bw()+th+colcale_cols_manual+stat_smooth(method='loess', color="black")
ggsave(file=here("figs","DPT_primary_PLA2G2A.pdf"),w=8, h=6)
ggsave(file=here("figs/jpeg","DPT_primary_PLA2G2A.jpeg"),w=8, h=6)
ggplot(plt_celltype, aes(pseudotime_diffusionmap, FABP1)) +geom_point(aes(colour = cluster_ID))+theme_bw()+th+colcale_cols_manual+stat_smooth(method='loess', color="black")
ggsave(file=here("figs","DPT_primary_FABP1.pdf"),w=8, h=6)
ggsave(file=here("figs/jpeg","DPT_primary_FABP1.jpeg"),w=8, h=6)

