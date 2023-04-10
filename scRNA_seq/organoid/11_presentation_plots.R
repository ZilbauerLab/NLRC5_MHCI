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
            #source(here("../../../codon/scRNAseq_codon/scripts","00_pretty_plots.R"))

##########################################
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

d10x.stim_norm<-readRDS(here("data","d10x_normalized.rds"))
identical(rownames(d10x.stim_norm@meta.data),rownames(d10x.stim@meta.data))
d10x.stim_norm <- AddMetaData(d10x.stim_norm, metadata = score_data)

#########
## MHCI
#########
FeaturePlot(d10x.stim_norm, features = "MHCI_score1",reduction = "umap", split.by ="orig.ident", min.cutoff = "q9", pt.size=1)
FeaturePlot(d10x.stim_norm, features = "crypt_villis_score1",reduction = "pca",  min.cutoff = "q9", pt.size=1)

umap_mat<-as.data.frame(Embeddings(object = d10x.stim_norm, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.stim_norm@meta.data
meta$cell<-rownames(meta)
rm(d10x.stim_norm)
gc()

plt<-merge(meta, umap_mat, by="cell")

plt$cluster_ID<-as.factor(plt$cluster_ID)
levels(plt$cluster_ID)<-c("Stem","Enterocyte","TA")

## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}



cluster_ID_plt_UMAP<-ggplot(plt, aes(UMAP_1,UMAP_2, color=cluster_ID))+
  geom_point(size=1.5)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  scale_color_manual(values=c("#6baed6","#238b45","#f16913"), name="Cell Type")
leg_umap_clust = get_leg(cluster_ID_plt_UMAP)

cluster_ID_plt_UMAP<- grid.arrange(cluster_ID_plt_UMAP+theme(legend.position = "none"),
                                   leg_umap_clust, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("figs/jpeg","aligned_cell_label_UMAP.jpeg"),cluster_ID_plt_UMAP, w=6, h=5)


stim_plt_UMAP<-ggplot(plt, aes(UMAP_1,UMAP_2, color=orig.ident))+geom_point(size=1.5)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  scale_color_manual(values=c("cornflowerblue","grey80","firebrick4"), name="Treatment")
leg_umap_stim = get_leg(stim_plt_UMAP)

stim_plt_UMAP<- grid.arrange(stim_plt_UMAP+theme(legend.position = "none"),
                             leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("figs/jpeg","aligned_stim_UMAP.jpeg"),stim_plt_UMAP, w=6, h=5)


MHC_plt_UMAP<-ggplot(plt, aes(UMAP_1,UMAP_2, color=MHCI_score1))+geom_point(size=1.5)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  labs(color = "Pathway\nActivation\nScore")
leg_umap_stim = get_leg(MHC_plt_UMAP)

MHC_plt_UMAP<- grid.arrange(MHC_plt_UMAP+theme(legend.position = "none"),
                             leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("figs/jpeg","aligned_MHC1_UMAP.jpeg"),MHC_plt_UMAP, w=6, h=5)



crypt_villis_plt_UMAP<-ggplot(plt, aes(UMAP_1,UMAP_2, color=crypt_villis_score1))+geom_point(size=1.5)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  labs(color = "Pathway\nActivation\nScore")
leg_umap_stim = get_leg(crypt_villis_plt_UMAP)

crypt_villis_plt_UMAP<- grid.arrange(crypt_villis_plt_UMAP+theme(legend.position = "none"),
                            leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("figs/jpeg","aligned_crypt_villis_UMAP.jpeg"),crypt_villis_plt_UMAP, w=6, h=5)





################################
## Aligned plot primary
################################
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
  name = 'MHCI_score'
)

d10x.primary <- AddModuleScore(
  object = d10x.primary,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)

score_data<-d10x.primary@meta.data[,c("MHCI_score1","crypt_villis_score1")]

## epithelial type labels
load(here("data","primary.epi.cells.RData"))

primary.epi.cells <- AddMetaData(primary.epi.cells, metadata = score_data)
#filter remaining immune (179 cells)
primary.epi.cells<-subset(primary.epi.cells, subset = cluster_ID != "B cell" & cluster_ID != "T cell"  & cluster_ID != "CD8 T cell")

d10x.primary_scores<-FetchData(object = d10x.primary, vars = c(MHCI))
d10x.primary_scores_epithelial<-d10x.primary_scores[which(rownames(d10x.primary_scores)%in%rownames(primary.epi.cells@meta.data)),]

d10x.primary_scores_epithelial$cell<-rownames(d10x.primary_scores_epithelial)
plt_counts<-melt(d10x.primary_scores_epithelial)

meta<-primary.epi.cells@meta.data
meta$cell<-rownames(meta)
plt_counts<-merge(plt_counts,meta[,c("cell","cluster_ID")] )



umap_mat<-as.data.frame(Embeddings(object = primary.epi.cells, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-primary.epi.cells@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")

plt_no_neo<-plt[which(plt$orig.ident!="Neonatal"),]
plt_no_neo<-plt_no_neo[which(plt_no_neo$cluster_ID!="Memory B cell"),]
plt_no_neo$cluster_ID<-as.character(plt_no_neo$cluster_ID)

plt_no_neo$cluster_ID<-as.factor(plt_no_neo$cluster_ID)
levels(plt_no_neo$cluster_ID)<-c("BEST4\nEnterocyte","Stem","Enterocyte","Enteroendocrine","Goblet","Paneth","Paneth (UC only)", "TA","Tuft","Paneth\n(UC only)")
cols_manual_less<-c( "#87d435ff",
                     "cornflowerblue",
                     "#238b43ff","#5e1e6fff","#aa0044ff",
                     "#d9c026ff", "#CFA218",
                     "#f16916ff","#a13da1ff")

## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}



cluster_ID_plt_UMAP<-ggplot(plt_no_neo, aes(UMAP_1,UMAP_2, color=cluster_ID))+
  geom_point(size=1.5)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  scale_color_manual(values=cols_manual_less, name="Cell Type")
leg_umap_clust = get_leg(cluster_ID_plt_UMAP)

cluster_ID_plt_UMAP<- grid.arrange(cluster_ID_plt_UMAP+theme(legend.position = "none"),
                                   leg_umap_clust, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("figs/jpeg","aligned_cell_label_UMAP_primary.jpeg"),cluster_ID_plt_UMAP, w=8, h=6.5)
ggsave(file=here("figs","aligned_cell_label_UMAP_primary.pdf"),cluster_ID_plt_UMAP, w=8, h=6.5)
#ggsave(file="../../../codon/scRNAseq_codon/figs/aligned_cell_label_UMAP.pdf",cluster_ID_plt_UMAP, w=8, h=6.5)


stim_plt_UMAP<-ggplot(plt_no_neo, aes(UMAP_1,UMAP_2, color=orig.ident))+geom_point(size=1.5)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  colscale_diagnosis
leg_umap_stim = get_leg(stim_plt_UMAP)

stim_plt_UMAP<- grid.arrange(stim_plt_UMAP+theme(legend.position = "none"),
                             leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("figs/jpeg","aligned_diagnosis_UMAP.jpeg"),stim_plt_UMAP, w=8, h=6.5)
ggsave(file=here("figs","aligned_diagnosis_UMAP.pdf"),stim_plt_UMAP, w=8, h=6.5)

MHC_plt_UMAP<-ggplot(plt_no_neo, aes(UMAP_1,UMAP_2, color=MHCI_score1))+geom_point(size=1.5)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  labs(color = "Pathway\nActivation\nScore")
leg_umap_stim = get_leg(MHC_plt_UMAP)

MHC_plt_UMAP<- grid.arrange(MHC_plt_UMAP+theme(legend.position = "none"),
                            leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("figs/jpeg","aligned_MHC1_UMAP_primary.jpeg"),MHC_plt_UMAP,  w=8, h=6.5)



crypt_villis_plt_UMAP<-ggplot(plt_no_neo, aes(UMAP_1,UMAP_2, color=crypt_villis_score1))+geom_point(size=1.5)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  labs(color = "Pathway\nActivation\nScore")
leg_umap_stim = get_leg(crypt_villis_plt_UMAP)

crypt_villis_plt_UMAP<- grid.arrange(crypt_villis_plt_UMAP+theme(legend.position = "none"),
                                     leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("figs/jpeg","aligned_crypt_villis_UMAP_primary.jpeg"),crypt_villis_plt_UMAP,  w=8, h=6.5)








####################
## Dot plots in each treatment
####################

# organoid markers
genes_df<-data.frame(gene=c("LGR5","ASCL2","FABP1","KRT19","HELLS","PCNA","MKI67","TOP2A", "FCGBP","SPINK4", "SOX4","CD24"),
                     cellType=c("Stem","Stem","Early E","Early E","TA1","TA1","TA2","TA2","Goblet","Goblet","Tuft","Tuft"))

rect<-ggplot()+geom_rect(aes(xmin=1:nrow(genes_df), xmax=1:nrow(genes_df)+1, ymin=1, ymax=2, 
                             fill=cellType),genes_df, color="black") + geom_hline(yintercept=-1, color="white")+
  theme_bw()+theme(legend.position = "none",
                   axis.title=element_blank(),
                   axis.text=element_blank(),
                   axis.ticks=element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_blank(),
                   panel.border = element_blank(),plot.margin = unit(c(0,5,0,2.4), "cm"))+
  geom_text(aes(x=seq(1,nrow(genes_df), 2)+1, y=0, label=genes_df$cellType[seq(1,nrow(genes_df), 2)]), size=4)+
  scale_fill_manual(values=c("#cf92beff","#364557","#9dbc26ff","#efe532ff","#F6A317","#89DFE2"))#





#'## Controls 
NT.cells.integrated<-readRDS(here("data","NT_normalized_integrated.rds"))
NT.cells.integrated

## add cell type labels from split analysis
load(here("output","celltype_labels_organoid.RData"))
cell_typelabel<-cell_typelabel[which(rownames(cell_typelabel)%in%rownames(NT.cells.integrated@meta.data)),]

identical(rownames(cell_typelabel),rownames(NT.cells.integrated@meta.data))
NT.cells.integrated <- AddMetaData(NT.cells.integrated, metadata = cell_typelabel)

marker_plt<-DotPlot(NT.cells.integrated, features = genes_df$gene, group.by='cluster_ID')+
  RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                          plot.margin = unit(c(1,0.5,0.1,0.1), "cm"))

ggsave(grid.arrange(marker_plt,rect,heights=c(0.9,0.1)),
       file=here("figs","NT_markerdotplot_celllabelled_SCT.pdf"), w=7,h=4)

ggsave(grid.arrange(marker_plt,rect,heights=c(0.9,0.1)),
       file=here("figs/jpeg","NT_markerdotplot_celllabelled_SCT.jpeg"), w=7,h=4)



#'## TNFa 
TNFa.cells.integrated<-readRDS(here("data","TNFa_normalized_integrated.rds"))
TNFa.cells.integrated

## add cell type labels from split analysis
load(here("output","celltype_labels_organoid.RData"))
cell_typelabel<-cell_typelabel[which(rownames(cell_typelabel)%in%rownames(TNFa.cells.integrated@meta.data)),]

identical(rownames(cell_typelabel),rownames(TNFa.cells.integrated@meta.data))
TNFa.cells.integrated <- AddMetaData(TNFa.cells.integrated, metadata = cell_typelabel)

marker_plt<-DotPlot(TNFa.cells.integrated, features = genes_df$gene, group.by='cluster_ID')+
  RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                          plot.margin = unit(c(1,0.5,0.1,0.1), "cm"))
# organoid markers noFCGBP
genes_dfnoFCGBP<-data.frame(gene=c("LGR5","ASCL2","FABP1","KRT19","HELLS","PCNA","MKI67","TOP2A", "SPINK4", "SOX4","CD24"),
                            cellType=c("Stem","Stem","Early E","Early E","TA1","TA1","TA2","TA2","Goblet","Tuft","Tuft"))

rectnoFCGBP<-ggplot()+geom_rect(aes(xmin=1:nrow(genes_dfnoFCGBP), xmax=1:nrow(genes_dfnoFCGBP)+1, ymin=1, ymax=2, 
                                    fill=cellType),genes_dfnoFCGBP, color="black") + geom_hline(yintercept=-1, color="white")+
  theme_bw()+theme(legend.position = "none",
                   axis.title=element_blank(),
                   axis.text=element_blank(),
                   axis.ticks=element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_blank(),
                   panel.border = element_blank(),plot.margin = unit(c(0,5,0,2.4), "cm"))+
  geom_text(aes(x=c(2,4,6,8,9.5,11), y=0, label=genes_dfnoFCGBP$cellType[seq(1,nrow(genes_dfnoFCGBP), 2)]), size=4)+
  scale_fill_manual(values=c("#cf92beff","#364557","#9dbc26ff","#efe532ff","#F6A317","#89DFE2"))#


ggsave(grid.arrange(marker_plt,rectnoFCGBP,heights=c(0.9,0.1)),
       file=here("figs","TNFa_markerdotplot_celllabelled_SCT.pdf"), w=7,h=4)

ggsave(grid.arrange(marker_plt,rectnoFCGBP,heights=c(0.9,0.1)),
       file=here("figs/jpeg","TNFa_markerdotplot_celllabelled_SCT.jpeg"), w=7,h=4)

#'## IFNg 
IFNg.cells.integrated<-readRDS(here("data","IFNg_normalized_integrated.rds"))
IFNg.cells.integrated

## add cell type labels from split analysis
load(here("output","celltype_labels_organoid.RData"))
cell_typelabel<-cell_typelabel[which(rownames(cell_typelabel)%in%rownames(IFNg.cells.integrated@meta.data)),]

identical(rownames(cell_typelabel),rownames(IFNg.cells.integrated@meta.data))
IFNg.cells.integrated <- AddMetaData(IFNg.cells.integrated, metadata = cell_typelabel)

marker_plt<-DotPlot(IFNg.cells.integrated, features = genes_df$gene, group.by='cluster_ID')+
  RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                          plot.margin = unit(c(1,0.5,0.1,0.1), "cm"))

ggsave(grid.arrange(marker_plt,rect,heights=c(0.9,0.1)),
       file=here("figs","IFNg_markerdotplot_celllabelled_SCT.pdf"), w=7,h=4)

ggsave(grid.arrange(marker_plt,rect,heights=c(0.9,0.1)),
       file=here("figs/jpeg","IFNg_markerdotplot_celllabelled_SCT.jpeg"), w=7,h=4)


print("## Organoid Dot plots combined")
####################
## Organoid Dot plots combined
####################

# organoid markers
genes_df<-data.frame(gene=c("LGR5","ASCL2","FABP1","KRT19","HELLS","PCNA","MKI67","TOP2A"),
                     cellType=c("Stem","Stem","Enterocyte","Enterocyte","TA1","TA1","TA2","TA2"))

rect<-ggplot()+geom_rect(aes(xmin=1:nrow(genes_df), xmax=1:nrow(genes_df)+1, ymin=1, ymax=2, 
                             fill=cellType),genes_df, color="black") + geom_hline(yintercept=-1, color="white")+
  theme_bw()+theme(legend.position = "none",
                   axis.title=element_blank(),
                   axis.text=element_blank(),
                   axis.ticks=element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_blank(),
                   panel.border = element_blank(),plot.margin = unit(c(0,5,0,2.4), "cm"))+
  geom_text(aes(x=seq(1,nrow(genes_df), 2)+1, y=0, label=genes_df$cellType[seq(1,nrow(genes_df), 2)]), size=4)+
  scale_fill_manual(values=c("#238b45","#6baed6","#d85e11","#f38742"))#f16913

## add cell type labels from split analysis
load(here("output","celltype_labels_organoid.RData"))
identical(rownames(cell_typelabel),rownames(d10x.stim@meta.data))
cell_typelabel$cluster_ID[which(cell_typelabel$cluster_ID=="entero_stem")]<-"enterocyte_stem"
d10x.stim <- AddMetaData(d10x.stim, metadata = cell_typelabel)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.stim <- NormalizeData(d10x.stim,scale.factor = 10000, normalization.method = "LogNormalize")

d10x.stim$cluster_ID<-as.factor(d10x.stim$cluster_ID)
levels(d10x.stim$cluster_ID)<-c("Stem", "Enterocyte","TA")

marker_plt<-DotPlot(d10x.stim, features = genes_df$gene, group.by='cluster_ID')+
  RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cell Type")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = unit(c(1,0.5,0.1,0.1), "cm"))+th_present

ggsave(grid.arrange(marker_plt,rect,heights=c(0.9,0.1)),
       file=here("figs","Allsamples_markerdotplot_celllabelled.pdf"), w=7,h=4)

ggsave(grid.arrange(marker_plt,rect,heights=c(0.9,0.1)),
       file=here("figs/jpeg","Allsamples_markerdotplot_celllabelled.jpeg"), w=7,h=4)



#' 
#' ####################
#' ## Primary Dot plots all cells
#' ####################
#' 
#' # primary markers
#' genes_df<-read.csv(here("data","single_cell_markers.csv"))
#' 
#' rect<-ggplot()+geom_rect(aes(xmin=1:nrow(genes_df), xmax=1:nrow(genes_df)+1, ymin=1, ymax=2, 
#'                              fill=Cell.Type),genes_df, color="black") + geom_hline(yintercept=-1, color="white")+
#'   theme_bw()+theme(legend.position = "none",
#'                    axis.title=element_blank(),
#'                    axis.text=element_blank(),
#'                    axis.ticks=element_blank(),
#'                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#'                    panel.background = element_blank(), axis.line = element_blank(),
#'                    panel.border = element_blank(),plot.margin = unit(c(0,5,0,2.4), "cm"))+
#'   geom_text(aes(x=seq(1,nrow(genes_df), 2)+1, y=0, label=genes_df$Cell.Type[seq(1,nrow(genes_df), 2)]), size=4)+
#'   scale_fill_manual(values=c("#238b45","#6baed6","#d85e11","#f38742"))#f16913
#' 
#' # compartment markers
#' rect_compartment<-ggplot()+geom_rect(aes(xmin=1:nrow(genes_df), xmax=1:nrow(genes_df)+1, ymin=1, ymax=2, 
#'                              fill=Compartment),genes_df, color="black") + geom_hline(yintercept=-1, color="white")+
#'   theme_bw()+theme(legend.position = "none",
#'                    axis.title=element_blank(),
#'                    axis.text=element_blank(),
#'                    axis.ticks=element_blank(),
#'                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#'                    panel.background = element_blank(), axis.line = element_blank(),
#'                    panel.border = element_blank(),plot.margin = unit(c(0,5,0,2.4), "cm"))+
#'   geom_text(aes(x=c(7,20,39), y=0, label=unique(genes_df$Compartment)), size=4)+
#'   scale_fill_manual(values=c("#238b45","#6baed6","#d85e11","#f38742"))#f16913
#' 
#' 
#' ## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
#' # but not normalized or scaled
#' d10x.primary<-readRDS("../../data/d10x_primary_raw_merged.rds")
#' 
#' ##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
#' # This is log(TP10K+1)
#' d10x.primary <- NormalizeData(d10x.primary,scale.factor = 10000, normalization.method = "LogNormalize")
#' 
#' ## add cell type labels from split analysis
#' #immune
#' load(here("../../output/","immune_iterative_label.Rdata"))
#' #stromal
#' load(here("../../output/","strom_celltype_label.Rdata"))
#' #epithelial
#' load(here("../../output/","epi_celltype_label.Rdata"))
#' cell_label<-rbind(epi_cell_labels, immune_cell_labels, stromal_cell_labels)
#' cell_label$cluster_ID<-as.character(cell_label$cluster_ID)
#' cell_label$cluster_ID[which(cell_label$cluster_ID=="Neonatal_cell")]<-"Neonatal Epithelial"
#' cell_label$index<-rownames(cell_label)
#' cell_label<-cell_label[match(colnames(d10x.primary), cell_label$index),]
#' identical(colnames(d10x.primary), cell_label$index)
#' 
#' d10x.primary <- AddMetaData(d10x.primary, metadata = cell_label)
#' 
#' table(d10x.primary@meta.data$cluster_ID)
#' 
#' d10x.primary@meta.data$cell_type<-as.factor(d10x.primary@meta.data$cluster_ID)
#' d10x.primary@meta.data$cell_type<-factor(d10x.primary@meta.data$cell_type, levels = c("Unidentified_only_UC",
#'                                                                                       "BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine",
#'                                                                                       "Goblet cell", "Paneth cell","Tuft",
#'                                                                                       "Activated B cell","activated DC" ,"Activated T","Arterial endothelial cell" ,"B cell","CD4 T cell",
#'                                                                                       "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
#'                                                                                       "IgA plasma cell","IgG plasma cell","Lymphatic endothelial cell","Macrophage" ,"mast cells","Memory B cell" ,
#'                                                                                       "Monocyte","pDC" ,"pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"TA","Tfh" ,
#'                                                                                       "Treg", "Venous endothelial cell","Glial cell","myofibroblast","Mast cell","unknown_filtered","Doublet","Neonatal CD4 T cell","Neonatal B cell","Neonatal Epithelial"))
#' 
#' #' ### dev cell labels
#' cols_manual<-c(  "pink",
#'                  "#3D0AF2","#6521C1","#67D364","#367C34",
#'                  "#B921C1","#308AC8",
#'                  "#C86E30","#C12134","#238443","#810f7c" ,"#02818a","#8c510a" ,"#78c679","#67a9cf",
#'                  "#3690c0","#8073ac","#b2abd2","#238443","#dd3497","#1d91c0","#addd8e","#014636",
#'                  "#253494","#081d58","#bf812d","#7a0177" ,"#ce1256","#006d2c" ,
#'                  "#a50f15","#542788" ,"#e31a1c","#e08214","#ef6548","#fdd49e" ,"#f46d43","#4575b4" ,
#'                  "#016c59", "#bf812d","#c51b7d","#fb9a99","#fb6a4a","lightgrey","black","#c6dbef","#c7e9c0","pink")
#' 
#' 
#' names(cols_manual) <- levels(d10x.primary@meta.data$cell_type)
#' fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
#' colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)
#' 
#' # ## remove doublets
#' d10x.primary<-subset(d10x.primary, subset = cluster_ID != "Doublet")
#' d10x.primary<-subset(d10x.primary, subset = cluster_ID != "Neonatal CD4 T cell")
#' d10x.primary<-subset(d10x.primary, subset = cluster_ID != "Neonatal Epithelial")
#' d10x.primary<-subset(d10x.primary, subset = cluster_ID != "Neonatal B cell")
#' 
#' 
#' genes_df<-genes_df[which(genes_df$Marker!="CD16"),]
#' 
#' marker_plt<-DotPlot(d10x.primary, features = unique(genes_df$Marker), group.by='cell_type')+
#'   RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cell Type")+
#'   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#'         plot.margin = unit(c(1,0.5,0.1,0.1), "cm"))+th_present
#' 
#' ggsave(grid.arrange(marker_plt,rect,heights=c(0.9,0.1)),
#'        file=here("figs","Allsamples_markerdotplot_celllabelled_primary.pdf"), w=7,h=4)
#' 
#' ggsave(grid.arrange(marker_plt,rect,heights=c(0.9,0.1)),
#'        file=here("figs/jpeg","Allsamples_markerdotplot_celllabelled_primary.jpeg"), w=7,h=4)
#' 
#' 

print("dot plot primary just epithelial")
################
## dot plot primary just epithelial
################
# organoid markers
genes_df<-data.frame(gene=c("LGR5","ASCL2","FABP1","KRT19","HELLS","PCNA","MKI67","TOP2A","FCGBP","SPINK4","SOX4","CD24"),
                     cellType=c("Stem","Stem","Enterocyte","Enterocyte","TA1","TA1","TA2","TA2","Goblet","Goblet","Tuft","Tuft"))

rect<-ggplot()+geom_rect(aes(xmin=1:nrow(genes_df), xmax=1:nrow(genes_df)+1, ymin=1, ymax=2, 
                             fill=cellType),genes_df, color="black") + geom_hline(yintercept=-1, color="white")+
  theme_bw()+theme(legend.position = "none",
                   axis.title=element_blank(),
                   axis.text=element_blank(),
                   axis.ticks=element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_blank(),
                   panel.border = element_blank(),plot.margin = unit(c(0,5,0,3.9), "cm"))+
  geom_text(aes(x=seq(1,nrow(genes_df), 2)+1, y=0, label=genes_df$cellType[seq(1,nrow(genes_df), 2)]), size=3.5)+
  scale_fill_manual(values=c("#238b43ff","#aa0044ff","cornflowerblue","#d85e11","#f38742","#a13da1ff"))#"#d85e11","#f38742"


## primary epithelial
## load raw primay data 
d10x.primary.raw<-readRDS(here("data","d10x_primary_raw_merged.rds"))

## pull comparment assignments from combined primary clustering/marker identification
load(here("output","cell_label_whole_clustering.RData"))

## add cell label information
cell_label$index<-rownames(cell_label)
cell_label<-cell_label[match(colnames(d10x.primary.raw), cell_label$index),]
identical(colnames(d10x.primary.raw), cell_label$index)
identical(rownames(d10x.primary.raw@meta.data), cell_label$index)

d10x.primary.raw<- AddMetaData(d10x.primary.raw, cell_label$seurat_clusters, col.name = "seurat_wholePrimary_clusters")
d10x.primary.raw<- AddMetaData(d10x.primary.raw, cell_label$elmentaite_annotation_V2, col.name = "elmentaite_annotation_V2")


## subset to just epithelial
primary.epi.cells<-subset(d10x.primary.raw, subset = seurat_wholePrimary_clusters %in% c(12,15,26))

primary.epi.cells_subset<-subset(primary.epi.cells, subset = individual %in% c(
  "T017","T019","T176","T189","T197","T202","T203","T024","T036","T44","T057",
  "T160","T161","T175","T182","T184","T180"))
rm(primary.epi.cells)

load(here("output","epi_celltype_label.Rdata"))

epi_cell_labels$index<-rownames(epi_cell_labels)
epi_cell_labels<-epi_cell_labels[match(colnames(primary.epi.cells_subset), epi_cell_labels$index),]
identical(colnames(primary.epi.cells_subset), epi_cell_labels$index)
epi_cell_labels$cluster_ID[which(epi_cell_labels$cluster_ID=="crypt")]<-"Stem"
epi_cell_labels$cluster_ID[which(epi_cell_labels$cluster_ID=="Unidentified_only_UC")]<-"Paneth (UC only)"


primary.epi.cells_subset <- AddMetaData(primary.epi.cells_subset, metadata = epi_cell_labels)
primary.epi.cells_subset<-subset(primary.epi.cells_subset, subset = cluster_ID != "Neonatal_cell")


##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
primary.epi.cells_subset <- NormalizeData(primary.epi.cells_subset,scale.factor = 10000, normalization.method = "LogNormalize")
head(primary.epi.cells_subset@meta.data)

marker_plt<-DotPlot(primary.epi.cells_subset, features = genes_df$gene, group.by='cluster_ID')+
  RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cell Type")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.margin = unit(c(1,0.5,0.1,0.1), "cm"))+th_present

grid.arrange(marker_plt,rect,heights=c(0.9,0.1))

ggsave(grid.arrange(marker_plt,rect,heights=c(0.9,0.1)),
       file=here("figs","Allsamples_markerdotplot_celllabelled_primary.pdf"), w=7,h=4)

ggsave(grid.arrange(marker_plt,rect,heights=c(0.9,0.1)),
       file=here("figs/jpeg","Allsamples_markerdotplot_celllabelled_primary.jpeg"), w=7,h=4)
#ggsave(grid.arrange(marker_plt,rect,heights=c(0.9,0.1)),file=here("../../../codon/scRNAseq_codon/figs","Allsamples_markerdotplot_celllabelled_primary.pdf"), w=7,h=4)



####################
## Score correlation plot  organoid
####################
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

d10x.stim_norm<-readRDS(here("data","d10x_normalized.rds"))
identical(rownames(d10x.stim_norm@meta.data),rownames(d10x.stim@meta.data))
d10x.stim_norm <- AddMetaData(d10x.stim_norm, metadata = score_data)

#### data frame for plotting
umap_mat<-as.data.frame(Embeddings(object = d10x.stim_norm, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-d10x.stim_norm@meta.data
meta$cell<-rownames(meta)
rm(d10x.stim_norm)
gc()

plt<-merge(meta, umap_mat, by="cell")

plt$cluster_ID<-as.factor(plt$cluster_ID)
levels(plt$cluster_ID)<-c("Stem","Enterocyte","TA")

scatter<-ggplot(plt, aes(crypt_villis_score1, MHCI_score1, fill=orig.ident, color=orig.ident))+geom_point(shape=21, color="black", size=1.5)+
  scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"))+scale_color_manual(values=c("cornflowerblue","grey80","#cb2626"))+
  theme_bw()+th_present+theme(legend.position = "none",plot.margin = unit(c(0,0,3.67,2.8), "cm"))+
  stat_smooth(method="lm",se=F)+xlab("Crypt-Villus Score")+ylab("MHC I Score")


violin_celltype<-ggplot(plt, aes(reorder(cluster_ID, crypt_villis_score1), crypt_villis_score1, fill=cluster_ID))+
  geom_violin()+scale_x_discrete(position = "bottom") +coord_flip()+
  scale_fill_manual(values=c("#6baed6","#238b45","#f16913"))+
  theme_bw()+th_present+ylab("")+xlab("")+theme(axis.title.x=element_blank(),
                                                axis.text.x=element_blank(),
                                                axis.ticks.x=element_blank(),
                                                legend.position = "none",
                                                plot.margin = unit(c(1.1,0,0,1), "cm"))

box_diagnosis<-ggplot(plt, aes(orig.ident, MHCI_score1, fill=orig.ident))+geom_boxplot()+
  scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"), name="Treatment")+xlab("Treatment")+
  theme_bw()+th_present+theme(axis.title.y =element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank(),
                              plot.margin = unit(c(0,0,3,0), "cm"),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mn_cryptvillus<-tapply(plt$crypt_villis_score1, plt$cluster_ID, mean)
cryptvillus_order<-names(mn_cryptvillus)[(order(mn_cryptvillus))]

plt_count<-melt(table(plt$orig.ident, plt$cluster_ID))
colnames(plt_count)<-c("orig.ident","cluster_ID","count")
plt_count$cluster_ID<-factor(plt_count$cluster_ID, cryptvillus_order)
placeholder<-ggplot(plt_count, aes(orig.ident, cluster_ID))+geom_tile(color="grey50", fill="white",size=0.4)+
  geom_text(aes(label=count),size=3.75,color="grey40")+ggtitle(" Cell Number")+
  theme_void()+theme(axis.title=element_blank(),
                     axis.text=element_blank(),
                     axis.ticks=element_blank(),
                     plot.margin = unit(c(0.4,2.88,0,0), "cm"),
                     plot.title = element_text(size = 12))
placeholder

grid_fig<-plot_grid(violin_celltype, placeholder, scatter, box_diagnosis, 
                    #align = 'hv',
                    ncol = 2, axis="lr", 
                    rel_heights = c(0.5,2),
                    rel_widths = c(3.75,1.25))

ggsave(file=here("figs", "MHCI_crpyt_summary_organoids_present.pdf"),grid_fig, w=10,h=10)
ggsave(file=here("figs/jpeg", "MHCI_crpyt_summary_organoids_present.jpeg"),grid_fig, w=10,h=10)





####################
## Score correlation plot primary
####################
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
  name = 'MHCI_score'
)

d10x.primary <- AddModuleScore(
  object = d10x.primary,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)

score_data<-d10x.primary@meta.data[,c("MHCI_score1","crypt_villis_score1")]

## epithelial type labels
load(here("data","primary.epi.cells.RData"))

primary.epi.cells <- AddMetaData(primary.epi.cells, metadata = score_data)
#filter remaining immune (179 cells)
primary.epi.cells<-subset(primary.epi.cells, subset = cluster_ID != "B cell" & cluster_ID != "T cell"  & cluster_ID != "CD8 T cell")

d10x.primary_scores<-FetchData(object = d10x.primary, vars = c(MHCI))
d10x.primary_scores_epithelial<-d10x.primary_scores[which(rownames(d10x.primary_scores)%in%rownames(primary.epi.cells@meta.data)),]

d10x.primary_scores_epithelial$cell<-rownames(d10x.primary_scores_epithelial)
plt_counts<-melt(d10x.primary_scores_epithelial)

meta<-primary.epi.cells@meta.data
meta$cell<-rownames(meta)
plt_counts<-merge(plt_counts,meta[,c("cell","cluster_ID")] )

cols_manual_less<-c( "#87d435ff",
                     "cornflowerblue",
                     "#238b43ff","#5e1e6fff","#aa0044ff",
                     "#d9c026ff","#CFA218",  
                     "#f16916ff","#a13da1ff")

umap_mat<-as.data.frame(Embeddings(object = primary.epi.cells, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-primary.epi.cells@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")

plt_no_neo<-plt[which(plt$orig.ident!="Neonatal"),]
plt_no_neo<-plt_no_neo[which(plt_no_neo$cluster_ID!="Memory B cell"),]
plt_no_neo$cluster_ID<-as.character(plt_no_neo$cluster_ID)

plt_no_neo$cluster_ID<-as.factor(plt_no_neo$cluster_ID)
levels(plt_no_neo$cluster_ID)<-c("BEST4 Enterocyte","Stem","Enterocyte","Enteroendocrine","Goblet","Paneth","Paneth (UC only)","TA","Tuft")

scatter<-ggplot(plt_no_neo, aes(crypt_villis_score1, MHCI_score1, fill=orig.ident, color=orig.ident))+geom_point(shape=21, color="black", size=1.5)+
  fillscale_diagnosis+theme_bw()+th_present+theme(legend.position = "none",plot.margin = unit(c(0,0,1.95,3.65), "cm"))+scale_color_manual(values=c("dodgerblue3","grey50","darkgoldenrod3"))+
  stat_smooth(method="lm",se=F)+xlab("Crypt-Villus Score")+ylab("MHC I Score")


violin_celltype<-ggplot(plt_no_neo, aes(reorder(cluster_ID, crypt_villis_score1), crypt_villis_score1, fill=cluster_ID))+
  geom_violin()+scale_x_discrete(position = "bottom") +coord_flip()+
  scale_fill_manual(values=cols_manual_less)+
  theme_bw()+th_present+ylab("")+xlab("")+theme(axis.title.x=element_blank(),
                                        axis.text.x=element_blank(),
                                        axis.ticks.x=element_blank(),
                                        legend.position = "none",
                                        plot.margin = unit(c(1,0,0,1), "cm"))

plt_no_neo$orig.ident<-as.character(plt_no_neo$orig.ident)
table(plt_no_neo$orig.ident)
box_diagnosis<-ggplot(plt_no_neo, aes(orig.ident, MHCI_score1, fill=orig.ident))+geom_boxplot()+xlab("Diagnosis")+
  fillscale_diagnosis+theme_bw()+th_present+theme(axis.title.y =element_blank(),
                                          axis.text.y=element_blank(),
                                          axis.ticks.y=element_blank(),
                                          plot.margin = unit(c(0,0.2,1,0), "cm"),
                                          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

mn_cryptvillus<-tapply(plt_no_neo$crypt_villis_score1, plt_no_neo$cluster_ID, mean)
cryptvillus_order<-names(mn_cryptvillus)[(order(mn_cryptvillus))]

plt_count<-melt(table(plt_no_neo$orig.ident, plt_no_neo$cluster_ID))
colnames(plt_count)<-c("orig.ident","cluster_ID","count")
plt_count$cluster_ID<-factor(plt_count$cluster_ID, cryptvillus_order)
placeholder<-ggplot(plt_count, aes(orig.ident, cluster_ID))+geom_tile(color="grey50", fill="white",size=0.4)+
  geom_text(aes(label=count),size=3.75,color="grey40")+ggtitle(" Cell Number")+
  theme_void()+theme(axis.title=element_blank(),
                     axis.text=element_blank(),
                     axis.ticks=element_blank(),
                     plot.margin = unit(c(0.4,3.1,0,0), "cm"),
                     plot.title = element_text(size = 12))
placeholder

grid_fig<-plot_grid(violin_celltype, placeholder, scatter, box_diagnosis, 
                    #align = 'hv',
                    ncol = 2, axis="lr", 
                    rel_heights = c(1,2),
                    rel_widths = c(3.75,1.25))

ggsave(file=here("figs", "MHCI_crpyt_summary_present.pdf"),grid_fig, w=10,h=10)
ggsave(file=here("figs/jpeg", "MHCI_crpyt_summary_present.jpeg"),grid_fig, w=10,h=10)

