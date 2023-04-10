
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

## testing factor
d10x.primary.cd.ctrl$cell_stim<-paste(d10x.primary.cd.ctrl$cluster_ID, d10x.primary.cd.ctrl$orig.ident, sep = "_")
Idents(d10x.primary.cd.ctrl) <- "cell_stim"

table(d10x.primary.cd.ctrl$cluster_ID, d10x.primary.cd.ctrl$orig.ident)



#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

cell_types<-unique(d10x.primary.cd.ctrl$cluster_ID)

contrasts_celltype_stim<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x.primary.cd.ctrl$cell_stim[grep(cell_types[x],d10x.primary.cd.ctrl$cell_stim)], repeats.allowed = FALSE)}))

contrasts_celltype_stim

nrow(contrasts_celltype_stim)

contrasts_celltype_stim[37,]<-c("enterocyte_CD","enterocyte_Ctrl")
contrasts_celltype_stim

d10x.primary.cd.ctrl
# this is 922,965 tests across all comparisons (24,945 genes, 37 comparisons)

diff_exp_all<-lapply(1:nrow(contrasts_celltype_stim), function(x){
  de<-FindMarkers(d10x.primary.cd.ctrl, ident.1 = contrasts_celltype_stim[x,1], ident.2 = contrasts_celltype_stim[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  print(paste(contrasts_celltype_stim[x,1],"vs", contrasts_celltype_stim[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_stim[x,1]
  de$cell.2<-contrasts_celltype_stim[x,2]
  de})


diff_exp_all<-do.call(rbind, diff_exp_all)

save(diff_exp_all, file=here("data","primary_diff_genes.RData"))
#load(file=here("../../../codon/scRNAseq_codon/data","primary_diff_genes.RData"))



## was NLRC5 even tested (or two lowly expressed?)
diff_exp_all[grep("NLRC5",diff_exp_all$gene),]
## What MHC I genes were tested
diff_exp_all[which(diff_exp_all$gene%in%MHCI),]


# add a significant threshold here for adjusted as not all in these lists are
diff_exp_all<-diff_exp_all[which(diff_exp_all$p_val_adj<0.005),]

table(diff_exp_all$gene)[which(table(diff_exp_all$gene)>10)]





plots <- VlnPlot(d10x.primary.cd.ctrl, features = c("CD74", "PSMB9", "IRF1"), split.by = "orig.ident", group.by = "cluster_ID",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

plots <- VlnPlot(d10x.primary.cd.ctrl, features = c("SERPINA1", "FABP1", "B2M"), split.by = "orig.ident", group.by = "cluster_ID",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)



diff_exp_all[grep("TAP1",diff_exp_all$gene),]
diff_exp_all[grep("NLRC5",diff_exp_all$gene),]
diff_exp_all[grep("PSMB8",diff_exp_all$gene),]
diff_exp_all[grep("B2M",diff_exp_all$gene),]


plots <- VlnPlot(d10x.primary.cd.ctrl, features = c("NLRC5","PSMB8","B2M","TAP1"), split.by = "orig.ident", group.by = "cluster_ID",
                 pt.size = 0, combine = FALSE)
wrap_plots(plots = plots, ncol = 1)

ggsave(file=here("figs","differential_MHC1_primary.pdf"), w=15,h=20)
ggsave(file=here("figs/jpeg","differential_MHC1_primary.jpeg"), w=15,h=20)



lapply(MHCI,function(x) diff_exp_all[grep(x,diff_exp_all$gene),])

          # ########
          # print("NLRC5 zero")
          # ########
          # meta<-d10x.primary.cd.ctrl@meta.data
          # meta$cell<-rownames(meta)
          # 
          # NLRC5<-FetchData(object = d10x.primary.cd.ctrl, vars = "NLRC5")
          # NLRC5$cell<-rownames(NLRC5)
          # NLRC5<-merge(NLRC5, meta, by="cell")
          # 
          # NLRC5$zero<-"value"
          # NLRC5$zero[which(NLRC5$NLRC5==0)]<-"zero"
          # 
          # tapply(NLRC5$cell, list(NLRC5$cluster_ID,NLRC5$orig.ident,NLRC5$zero), function(x) length(unique(x)))

#############
#save expression for plotting
#############
meta<-d10x.primary.cd.ctrl@meta.data
meta$cell<-rownames(meta)

MHCI_exp<-FetchData(object = d10x.primary.cd.ctrl, vars = MHCI)
MHCI_exp$cell<-rownames(MHCI_exp)
MHCI_exp<-merge(MHCI_exp, meta, by="cell")
save(MHCI_exp, file=here("data","MHCI_primary_expression.RData"))
#load(here("../../../codon/scRNAseq_codon/data","MHCI_primary_expression.RData"))

MHCI_exp_epi<-MHCI_exp[which(MHCI_exp$cluster_ID%in%c("crypt","TA","BEST4 enterocyte","enterocyte","enteroendocrine","Goblet cell","Paneth")),]

MHCI_exp_epi<-melt(MHCI_exp_epi[,c( "cell","orig.ident","HLA-F","HLA-G","HLA-A","HLA-E","HLA-C","HLA-B","TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5","cluster_ID")], id=c("cluster_ID","cell","orig.ident"))

MHCI_exp_epi$cluster_ID<-as.factor(MHCI_exp_epi$cluster_ID)
levels(MHCI_exp_epi$cluster_ID)<-c("BEST4\nEnterocyte","Stem","Enterocyte","Entero-\nendocrine","Goblet","Paneth","TA")
MHCI_exp_epi$orig.ident<-as.factor(MHCI_exp_epi$orig.ident)
levels(MHCI_exp_epi$orig.ident)<-c("CD","Control")

diff_exp_all_MHC<-diff_exp_all[which(diff_exp_all$gene%in%MHCI),]
diff_exp_all_MHC$celltype<-sapply(1:nrow(diff_exp_all_MHC), function(x) strsplit(diff_exp_all_MHC$cell.1[x],"_")[[1]][1])
diff_exp_all_MHC_epithelial<-diff_exp_all_MHC[diff_exp_all_MHC$celltype%in%c("crypt","TA","BEST4 enterocyte","enterocyte","enteroendocrine","Goblet cell","Paneth" ),]

highlight<-diff_exp_all_MHC_epithelial[,c("gene","celltype")]
colnames(highlight)<-c("variable","cluster_ID")

highlight$cluster_ID<-as.factor(highlight$cluster_ID)
levels(highlight$cluster_ID)<-c("Stem","Enterocyte","TA")

        # ggplot()+geom_violin(aes(orig.ident, value, fill=orig.ident),MHCI_exp_epi)+facet_grid(variable~cluster_ID)+theme_bw()+th_present+
        #   scale_fill_manual(values=c("dodgerblue3","lightgrey"))+theme(legend.position = "none")+ylab("Expression Level")+xlab("Diagnosis")+
        #   geom_rect(data = highlight, 
        #           fill = NA, colour = "yellow", xmin = -Inf,xmax = Inf,
        #           ymin = -Inf,ymax = Inf, size=1.75)


ggplot()+geom_violin(aes(orig.ident, value, fill=orig.ident),MHCI_exp_epi, fill="grey", color="grey")+
  geom_boxplot(aes(orig.ident, value, fill=orig.ident),MHCI_exp_epi, width=0.25, outlier.shape=NA)+
  facet_grid(variable~cluster_ID, scales="free_y")+theme_bw()+th_present+
  scale_fill_manual(values=c("dodgerblue3","lightgrey"))+
  theme(legend.position = "none")+ylab("Expression Level")+xlab("Diagnosis")+
  geom_rect(data = highlight, 
            fill = NA, colour = "yellow", xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, size=1.75)
ggsave(file=here("figs","primary_MHC_singlecell.pdf"), w=8, h=10)
ggsave(file=here("figs/jpeg","primary_MHC_singlecell.jpeg"), w=8, h=10)










# ##############
# ## Plot genes
# ##############
# load(file=here("../../../codon/scRNAseq_codon/data","primary_diff_genes.RData"))
# diff_exp_all_MHC<-diff_exp_all[which(diff_exp_all$gene%in%MHCI),]
# 
# diff_exp_all_MHC$celltype<-sapply(1:nrow(diff_exp_all_MHC), function(x) strsplit(diff_exp_all_MHC$cell.1[x],"_")[[1]][1])
# diff_exp_all_MHC$significant<-sapply(1:nrow(diff_exp_all_MHC), function(x) if(diff_exp_all_MHC$p_val_adj[x]<0.005){"Significant"}else{"Not Significant"})
# 
# diff_exp_all_MHC_epithelial<-diff_exp_all_MHC[diff_exp_all_MHC$celltype%in%c("crypt","TA","BEST4 enterocyte","enterocyte","enteroendocrine","Goblet cell","Paneth" ),]
# 
# heat_plot_MHCI<-ggplot(diff_exp_all_MHC_epithelial, aes(celltype, gene, fill=avg_log2FC, size=significant))+geom_point()
# heat_plot_MHCI
# 
# 
# diff_exp_all_other<-diff_exp_all[which(diff_exp_all$gene%in%c("CIITA","HLA-DRA","PAX6")),]
# 
# diff_exp_all_other_INFgNT<-diff_exp_all_other[grep("IFNg",diff_exp_all_other$cell.1),]
# diff_exp_all_other_INFgNT<-diff_exp_all_other_INFgNT[grep("NT",diff_exp_all_other_INFgNT$cell.2),]
# 
# diff_exp_all_other_INFgNT$celltype<-sapply(1:nrow(diff_exp_all_other_INFgNT), function(x) strsplit(diff_exp_all_other_INFgNT$cell.1[x],"_")[[1]][1])
# diff_exp_all_other_INFgNT$significant<-sapply(1:nrow(diff_exp_all_other_INFgNT), function(x) if(diff_exp_all_other_INFgNT$p_val_adj[x]<0.005){"Significant"}else{"Not Significant"})
# 
# heat_plot_other<-ggplot(diff_exp_all_other_INFgNT, aes(celltype, gene, fill=avg_log2FC, size=significant))+geom_point()
# heat_plot_other
# 



#' ### venn
#' crypt_NT_IFNg<-diff_exp_all[which(diff_exp_all$cell.1=="crypt_IFNg" & diff_exp_all$cell.2=="crypt_NT"),]
#' TA_NT_IFNg<-diff_exp_all[which(diff_exp_all$cell.1=="TA_IFNg" & diff_exp_all$cell.2=="TA_NT"),]
#' entero_NT_IFNg<-diff_exp_all[which(diff_exp_all$cell.1=="enterocyte_IFNg" & diff_exp_all$cell.2=="enterocyte_NT"),]
#' enteromt_NT_IFNg<-diff_exp_all[which(diff_exp_all$cell.1=="enterocyte_mt_IFNg" & diff_exp_all$cell.2=="enterocyte_mt_NT"),]
#'
#'
#' intersect(crypt_NT_IFNg$gene, TA_NT_IFNg$gene)
#' intersect(crypt_NT_IFNg$gene, entero_NT_IFNg$gene)
#' intersect(crypt_NT_IFNg$gene, enteromt_NT_IFNg$gene)
#'
#' intersect(enteromt_NT_IFNg$gene, entero_NT_IFNg$gene)
#'
#' diff_inallcell<-intersect(intersect(TA_NT_IFNg$gene, entero_NT_IFNg$gene),intersect(crypt_NT_IFNg$gene, enteromt_NT_IFNg$gene))
#' intersect(MHCI,diff_inallcell)
#' 
#' 
#' crypt_specific<-crypt_NT_IFNg$gene[which(!(crypt_NT_IFNg$gene%in%entero_NT_IFNg$gene | crypt_NT_IFNg$gene%in%TA_NT_IFNg$gene | crypt_NT_IFNg$gene%in%enteromt_NT_IFNg$gene))]
#' entero_specific<-entero_NT_IFNg$gene[which(!(entero_NT_IFNg$gene%in%crypt_NT_IFNg$gene | entero_NT_IFNg$gene%in%TA_NT_IFNg$gene | entero_NT_IFNg$gene%in%enteromt_NT_IFNg$gene))]
#' enteromt_specific<-enteromt_NT_IFNg$gene[which(!(enteromt_NT_IFNg$gene%in%crypt_NT_IFNg$gene | enteromt_NT_IFNg$gene%in%TA_NT_IFNg$gene | enteromt_NT_IFNg$gene%in%entero_NT_IFNg$gene))]
#' TA_specific<-TA_NT_IFNg$gene[which(!(TA_NT_IFNg$gene%in%entero_NT_IFNg$gene | TA_NT_IFNg$gene%in%crypt_NT_IFNg$gene | TA_NT_IFNg$gene%in%enteromt_NT_IFNg$gene))]
#' 
#' entero_general<-entero_NT_IFNg$gene[which(!(entero_NT_IFNg$gene%in%crypt_NT_IFNg$gene | entero_NT_IFNg$gene%in%TA_NT_IFNg$gene))]
#' 
#' intersect(MHCI,crypt_specific)
#' intersect(MHCI,entero_specific)
#' intersect(MHCI,enteromt_specific)
#' intersect(MHCI,TA_specific)
#' 
#' intersect(MHCI,entero_general)
#' 
#' 
#' all_cell_genes_plt <- VlnPlot(d10x.stim, features = c("TAP1","PSMB9","B2M","IRF1"), split.by = "orig.ident", group.by = "cluster_ID", 
#'                  pt.size = 0, combine = FALSE)+scale_color_manual(values=c("cornflowerblue","grey80","firebrick4"))
#' wrap_plots(plots = all_cell_genes_plt, ncol = 1)
#' 
#' 
#' entero_genes_plt <- VlnPlot(d10x.stim, features = c("NLRC5","CIITA","MUC1","TRIM15"), split.by = "orig.ident", group.by = "cluster_ID", 
#'                               pt.size = 0, combine = FALSE)
#' wrap_plots(plots = entero_genes_plt, ncol = 1)
#' 
#' crypt_genes_plt <- VlnPlot(d10x.stim, features = c("MT-CO1","DUT","PSMB1","PSMB2"), split.by = "orig.ident", group.by = "cluster_ID", 
#'                             pt.size = 0, combine = FALSE)
#' wrap_plots(plots = crypt_genes_plt, ncol = 1)
#' 
#' TA_genes_plt <- VlnPlot(d10x.stim, features = c("HLA-DMB","ERAP2","SCML1","PCLAF"), split.by = "orig.ident", group.by = "cluster_ID", 
#'                            pt.size = 0, combine = FALSE)
#' wrap_plots(plots = TA_genes_plt, ncol = 1)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' ### gene plots
#'         #FetchData(object = d10x.stim, vars = c("UMAP_1", "UMAP_2", "MS4A1"))
#' 
#' 
#' 
#' 
#' ## normalized UMAP scores
#' d10x.stim_norm<-readRDS("../../data/d10x_normalized.rds")
#' 
#' umap_mat<-as.data.frame(Embeddings(object = d10x.stim_norm, reduction = "umap"))#
#' umap_mat$cell<-rownames(umap_mat)
#' 
#' meta<-d10x.stim_norm@meta.data
#' meta$cell<-rownames(meta)
#' rm(d10x.stim_norm)
#' gc()
#' # 
#' # ## raw expression values
#' # d10x.stim.exp<-as.data.frame(d10x.stim[["RNA"]]@data)
#' # rm(d10x.stim)
#' # 
#' # d10x.stim.exp.GOI<-d10x.stim.exp[unique(diff_exp_all$gene),]
#' # rm(d10x.stim.exp)
#' # gc()
#' # 
#' # d10x.stim.exp.GOI$gene<-rownames(d10x.stim.exp.GOI)
#' # d10x.stim.exp.GOI<-melt(d10x.stim.exp.GOI)#
#' 
#' 
#' ## merge to plot
#' plt<-merge(meta, umap_mat, by="cell")
#' cell_typelabel$cell<-rownames(cell_typelabel)
#' plt<-merge(plt, cell_typelabel, by="cell")
#' 
#' 
#' 
#' 
#' 
#' 
#' ## grab legened from plot
#' get_leg = function(a.gplot){
#'   tmp <- ggplot_gtable(ggplot_build(a.gplot))
#'   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#'   legend <- tmp$grobs[[leg]]
#'   legend
#' }
#' 
#' 
#' cluster_ID_plt_UMAP<-ggplot(plt, aes(UMAP_1,UMAP_2, color=cluster_ID))+
#'   geom_point(size=1.5)+
#'   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#'   scale_color_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))
#' leg_umap_clust = get_leg(cluster_ID_plt_UMAP)
#' 
#' cluster_ID_plt_UMAP<- grid.arrange(cluster_ID_plt_UMAP+theme(legend.position = "none"),
#'                               leg_umap_clust, ncol=2, widths=c(0.8,0.2))
#' ggsave(file="../../figs/scRNAseq/jpeg/cell_label_UMAP.jpeg",cluster_ID_plt_UMAP, w=6, h=5)
#' 
#' 
#' stim_plt_UMAP<-ggplot(plt, aes(UMAP_1,UMAP_2, color=orig.ident.x))+geom_point(size=1.5)+
#'   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#'   scale_color_manual(values=c("cornflowerblue","grey80","firebrick4"), name="Treatment")
#' leg_umap_stim = get_leg(stim_plt_UMAP)
#' 
#' stim_plt_UMAP<- grid.arrange(stim_plt_UMAP+theme(legend.position = "none"),
#'                              leg_umap_stim, ncol=2, widths=c(0.8,0.2))
#' 
#' 
#' 
#' 
#' 
#' 
#' ## Multiple plotting functions
#' gene_plot_UMAP<-function(gene_name){
#'   plt_gene<-FetchData(object = d10x.stim, vars = gene_name)
#'   plt_gene$cell<-rownames(plt_gene)
#'   plt_gene<-merge(plt_gene, plt, by="cell")
#'   plt_gene<-plt_gene[order(plt_gene[,2]),]
#' 
#'   ## UMAP to match others
#'   gene_plt_UMAP<-ggplot(plt_gene, aes(UMAP_1,UMAP_2, color=plt_gene[,2]))+geom_point(size=1.5)+
#'     theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
#'     scale_colour_gradient2( low = "#2b8cbe",#blue
#'                             mid = "#f0f0f0", # grey
#'                             high = "#eb8423", #oragne
#'                             midpoint = max(plt_gene[,gene_name])/2,
#'                             na.value = "#226d94",
#'                             name=paste(gene_name,"\nCount",sep=""))
#'   leg_umap_gene = get_leg(gene_plt_UMAP)
#'   
#'   grid.arrange(gene_plt_UMAP+theme(legend.position = "none"), leg_umap_gene, ncol=2, widths=c(0.8,0.2))
#'   }
#' 
#' 
#' 
#' gene_plot_UMAP_facet<-function(gene_name){
#'   plt_gene<-FetchData(object = d10x.stim, vars = gene_name)
#'   plt_gene$cell<-rownames(plt_gene)
#'   plt_gene<-merge(plt_gene, plt, by="cell")
#'   plt_gene<-plt_gene[order(plt_gene[,2]),]
#'   ggplot(plt_gene, aes(UMAP_1,UMAP_2, color=plt_gene[,2]))+
#'          geom_point(size=1.5)+
#'          scale_colour_gradient2( low = "#2b8cbe",#blue
#'                                  mid = "#f0f0f0", # grey
#'                                  high = "#eb8423", #oragne
#'                                  midpoint = max(plt_gene[,gene_name])/2,
#'                                  na.value = "#226d94",
#'                                  name="Count")+theme_classic()+
#'          ggtitle(gene_name)+th+theme(plot.title = element_text(size = 16, face = "bold"))+
#'          facet_grid(cluster_ID~orig.ident.x)
#'   }
#' 
#' gene_plot_box<-function(gene_name){
#'   plt_gene<-FetchData(object = d10x.stim, vars = gene_name)
#'   plt_gene$cell<-rownames(plt_gene)
#'   plt_gene<-merge(plt_gene, plt, by="cell")
#'   plt_gene<-plt_gene[order(plt_gene[,2]),]
#'   
#'   ## seperated boxplot
#'   plt_speprate_cell<-function(cell){
#'     plt_gene_cell<-plt_gene[which(plt_gene$cluster_ID==cell),]
#'     sig<-diff_exp_all[intersect(grep(gene_name,diff_exp_all$gene), intersect(grep(cell,diff_exp_all$cell.1), grep(cell,diff_exp_all$cell.2))),]
#'     
#'     if(nrow(sig)>0){
#'       if(cell=="enterocyte"){
#'         rm_rows<-unique(c(which(grepl("mt",sig$cell.1)), which(grepl("mt",sig$cell.2))))
#'         if(length(rm_rows)>0){sig<-sig[-rm_rows,]}}
#'       if(cell=="enterocyte_mt"){
#'         comp<-lapply(1:nrow(sig), function(x){c(strsplit(sig[x,"cell.1"],"_")[[1]][3], strsplit(sig[x,"cell.2"],"_")[[1]][3])})}else{
#'           comp<-lapply(1:nrow(sig), function(x){c(strsplit(sig[x,"cell.1"],"_")[[1]][2], strsplit(sig[x,"cell.2"],"_")[[1]][2])})
#'           }
#'       
#'       ggplot(plt_gene_cell, aes(orig.ident.x, plt_gene_cell[,2], fill=orig.ident.x))+geom_violin()+
#'         scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"), guide=F)+
#'         geom_signif(comparisons =comp,
#'                     step_increase = 0.1,map_signif_level = T, test=NULL, annotations=signif(sig$p_val_adj,3))+
#'         theme_bw()+th+xlab("")+ylab(paste(gene_name, "Count"))+
#'         ggtitle(cell)+ylim(0,round(max(plt_gene[,2]))+(round(max(plt_gene[,2]))/2.5))
#'       }else{
#'           ggplot(plt_gene_cell, aes(orig.ident.x, plt_gene_cell[,2], fill=orig.ident.x))+geom_violin()+
#'             scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"), guide=F)+
#'             theme_bw()+th+xlab("")+ylab(paste(gene_name, "Count"))+
#'             ggtitle(cell)+ylim(0,round(max(plt_gene[,2]))+(round(max(plt_gene[,2]))/2.5))
#'         }
#'   }
#'   
#'   leg = get_leg(ggplot(plt_gene, aes(orig.ident.x, plt_gene[,2], fill=orig.ident.x))+geom_violin()+
#'                   scale_fill_manual(values=c("cornflowerblue","grey80","firebrick4"), name="Treatment"))
#'   
#'   grid.arrange(arrangeGrob(plt_speprate_cell("crypt"),
#'                                             plt_speprate_cell("TA"),
#'                                             plt_speprate_cell("enterocyte"),
#'                                             plt_speprate_cell("enterocyte_mt"),  ncol=2),
#'                                 leg, ncol=2, widths=c(0.8,0.2))
#'   }
#' 
#' 
#' 
#' # look at one gene quick
#' gene<-"NLRC5"
#' UMAP_gene<-gene_plot_UMAP(gene)
#' umaps_gene<-grid.arrange(stim_plt_UMAP, cluster_ID_plt_UMAP,UMAP_gene,ncol=2)
#' ggsave(file=paste("../../figs/scRNAseq/gene_plots/jpeg/",gene,"_UMAP.jpeg",sep=""),umaps_gene, w=11, h=10)
#' 
#' gene_plot_UMAP_facet(gene)
#' ggsave(file=paste("../../figs/scRNAseq/gene_plots/jpeg/",gene,"_UMAP_facet.jpeg",sep=""), w=10, h=9)
#' 
#' box_gene<-gene_plot_box(gene)
#' ggsave(file=paste("../../figs/scRNAseq/gene_plots/jpeg/",gene,"_box.jpeg",sep=""),box_gene, w=6, h=7)
#' 
#' 
#' #plot a bunch
#' lapply(c("NLRC5","TAP1","TAP2","B2M","CIITA","SERPINA1","IRF1",'HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"PSMB9","PSMB8","MR1","CD1D"), function(gene){
#'     UMAP_gene<-gene_plot_UMAP(gene)
#'     umaps_gene<-grid.arrange(stim_plt_UMAP, cluster_ID_plt_UMAP,UMAP_gene,ncol=2)
#'     ggsave(file=paste("../../figs/scRNAseq/gene_plots/jpeg/",gene,"_UMAP.jpeg",sep=""),umaps_gene,  w=11, h=10)
#'     ggsave(file=paste("../../figs/scRNAseq/gene_plots/",gene,"_UMAP.pdf",sep=""),umaps_gene, w=11, h=10)
#'     
#'     
#'     gene_plot_UMAP_facet(gene)
#'     ggsave(file=paste("../../figs/scRNAseq/gene_plots/jpeg/",gene,"_UMAP_facet.jpeg",sep=""), w=10, h=9)
#'     ggsave(file=paste("../../figs/scRNAseq/gene_plots/",gene,"_UMAP_facet.pdf",sep=""), w=10, h=9)
#'     
#'     
#'     box_gene<-gene_plot_box(gene)
#'     ggsave(file=paste("../../figs/scRNAseq/gene_plots/jpeg/",gene,"_box.jpeg",sep=""),box_gene, w=6, h=7)
#'     ggsave(file=paste("../../figs/scRNAseq/gene_plots/",gene,"_box.pdf",sep=""),box_gene, w=6, h=7)})
#' 
#' 
#' 
#' 
#' 
#' 
#' ##################
#' #'## Genes differential Entero NT and entero IFNg
#' ##################
#' entero_NT_IFNg<-diff_exp_all[which(diff_exp_all$cell.1=="enterocyte_IFNg" & diff_exp_all$cell.2=="enterocyte_NT"),]
#' box_gene<-gene_plot_box("FABP1")
#' box_gene<-gene_plot_box("MUC1")
#' 
#' #' unique different between entero
#' 
#' intersect(crypt_NT_IFNg$gene, TA_NT_IFNg$gene)
#' intersect(crypt_NT_IFNg$gene, entero_NT_IFNg$gene)
#' intersect(crypt_NT_IFNg$gene, enteromt_NT_IFNg$gene)
#' 
#' intersect(enteromt_NT_IFNg$gene, entero_NT_IFNg$gene)
#' 
#' diff_inallcell<-intersect(intersect(TA_NT_IFNg$gene, entero_NT_IFNg$gene),intersect(crypt_NT_IFNg$gene, enteromt_NT_IFNg$gene))
#' intersect(MHCI,diff_inallcell)
#' diff_inanycell<-unique(c(TA_NT_IFNg$gene, entero_NT_IFNg$gene,crypt_NT_IFNg$gene, enteromt_NT_IFNg$gene))
#' write.table(diff_inanycell, file=here("../../output/IFNgdiff_genes_any_cell_type.txt"), quote=F, row.names = F, col.names = F)
#' 
#' 
#' entero_NT_IFNg$gene[which(!(entero_NT_IFNg$gene %in% c(TA_NT_IFNg$gene, crypt_NT_IFNg$gene, enteromt_NT_IFNg$gene)))]
#' 
#' 
#' write.table(entero_NT_IFNg$gene, file=here("../../output/entero_IFNgdiff_genes.txt"), quote=F, row.names = F, col.names = F)
#' write.table(entero_NT_IFNg$gene[which(!(entero_NT_IFNg$gene %in% c(TA_NT_IFNg$gene, crypt_NT_IFNg$gene, enteromt_NT_IFNg$gene)))], file=here("../../output/entero_specific_IFNgdiff_genes.txt"), quote=F, row.names = F, col.names = F)
#' # run in ermineJ. Interferion stuff comes up with general list and no specific hormaone signalling comes up in entero specific list
#' 
#' # stress genes
#' gene_plot_box("DNAJB1")
#' gene_plot_box("HSPA1A")
#' diff_exp_all[which(diff_exp_all$gene=="DNAJB1"),]
#' diff_exp_all[which(diff_exp_all$gene=="HSPA1A"),]
#' 
#' 
#' 
#' 
#' 
#' ###############
#' ## Gene correlation
#' ##############
#' NLRC5<-FetchData(object = d10x.stim, vars = "NLRC5")
#' NLRC5$cell<-rownames(NLRC5)
#' NLRC5<-merge(NLRC5, plt, by="cell")
#' NLRC5<-NLRC5[order(NLRC5[,2]),]
#' 
#' gene_correlation<-function(gene_name){
#'   plt_gene<-FetchData(object = d10x.stim, vars = c("NLRC5",gene_name))
#'   plt_gene$cell<-rownames(plt_gene)
#'   plt_gene<-merge(plt_gene, plt, by="cell")
#'   
#'   ## ignore 0s for correlation?
#'   plt_gene$zero<-"value"
#'   plt_gene$zero[which(plt_gene$NLRC5==0 | plt_gene[,3]==0)]<-"zero"
#'   
#' split_cor<-split(plt_gene, list(plt_gene$cluster_ID, plt_gene$orig.ident.y))
#' 
#' do.call(rbind,lapply(split_cor, function(df){
#'   if(nrow(df)==0){}else{
#'     df<-df[which(df$zero=="value"),]
#'     data.frame(stim=unique(df$orig.ident.y), cell=unique(df$cluster_ID), gene=gene_name, correlation=cor(df$NLRC5, df[,3]))}
#' }))
#' 
#' ggplot(plt_gene, aes(NLRC5, plt_gene[,3],  color=zero))+geom_point()+ylab(gene_name)+
#'   facet_grid(cluster_ID~orig.ident.x)+theme_bw()+th+scale_color_manual(values=c("cornflowerblue","lightgrey"))
#' }
#' 
#' 
#' gene_correlation_value<-function(gene_name){
#'   plt_gene<-FetchData(object = d10x.stim, vars = c("NLRC5",gene_name))
#'   plt_gene$cell<-rownames(plt_gene)
#'   plt_gene<-merge(plt_gene, plt, by="cell")
#'   
#'   ## ignore 0s for correlation?
#'   plt_gene$zero<-"value"
#'   plt_gene$zero[which(plt_gene$NLRC5==0 | plt_gene[,3]==0)]<-"zero"
#'   paste(gene_name,":",cor(plt_gene$NLRC5,plt_gene[,3]))
#'   }
#'   
#' lapply(1:length(MHCI), function(x) gene_correlation_value(MHCI[x]))
#' gene_correlation("TAP1")
#' ggsave(file="../../figs/scRNAseq/gene_plots/jpeg/TAP1_correlation_w_NLRC5.jpeg", w=7, h=6)
#' ggsave(file="../../figs/scRNAseq/gene_plots/TAP1_correlation_w_NLRC5.pdf", w=7, h=6)
#' gene_correlation("IRF1")
#' ggsave(file="../../figs/scRNAseq/gene_plots/jpeg/IRF1_correlation_w_NLRC5.jpeg", w=7, h=6)
#' ggsave(file="../../figs/scRNAseq/gene_plots/IRF1_correlation_w_NLRC5.pdf", w=7, h=6)
#' 
#' 
#' 
#' 
#' 
#' ##################
#' # Heat map
#' ##################
#' genes<-unique(c(MHCI,diff_inallcell[1:5],crypt_specific[1:5],entero_specific[1:5],enteromt_specific[1:5],TA_specific[1:5]))
#' 
#' d10x.stim <- ScaleData(d10x.stim) 
#' DoHeatmap(subset(d10x.stim, downsample = 100), features = genes, size = 3)
#' 
#' 
#' DoHeatmap(d10x.stim, features=c(MHCI,diff_inallcell))
#' 
#' DoHeatmap(d10x.stim, features=c(MHCI,entero_NT_IFNg$gene[1:100]))
#' 
#' 
#' 
#' # 
#' # plt_gene<-FetchData(object = d10x.stim, vars = genes)
#' # plt_gene<-plt_gene[sample(1:nrow(plt_gene),5000),]
#' # plt_gene$cell<-rownames(plt_gene)
#' # plt_gene<-melt(plt_gene)#
#' # 
#' # plt_gene<-merge(plt_gene, plt, by="cell")
#' # 
#' # plt_gene<-plt_gene[order(plt_gene$orig.ident.x),]
#' # 
#' # 
#' # plt_gene_crypt<-plt_gene[which(plt_gene$cluster_ID=="crypt"),]
#' # 
#' # labels<-plt_gene_crypt[,c("cell","orig.ident.x")]
#' # labels<-labels[!duplicated(labels),]
#' # 
#' # ggplot(plt_gene_crypt, aes(cell, variable, fill= value)) +   geom_tile() + 
#' #   theme_classic()+ 
#' #   scale_x_discrete(labels=labels$orig.ident.x)#+
#' #   # theme(axis.title.x=element_blank(),
#' #   #       axis.text.x=element_blank(),
#' #   #       axis.ticks.x=element_blank())
#' 
#' 
#' 
#' ############
#' ## MHC I diff summary
#' ############
#' IFNg_NT_diff<-diff_exp_all[-unique(c(grep("TNFa", diff_exp_all$cell.1),grep("TNFa", diff_exp_all$cell.2))),]
#' IFNg_NT_diff_MHCI<-IFNg_NT_diff[which(IFNg_NT_diff$gene%in%MHCI),]
#' IFNg_NT_diff_MHCI$cell<-gsub("_IFNg|_NT","",IFNg_NT_diff_MHCI$cell.1)
#' 
#' ggplot(IFNg_NT_diff_MHCI, aes(gene, cell))+geom_tile(aes(fill=avg_logFC),color = "black",size=0.5)+
#'   geom_text(aes(label=round(avg_logFC,2)))+ 
#'   theme_gray(8)+
#'   theme(axis.text = element_text(size =12, color="black"),
#'         axis.title = element_text(size =15),
#'         legend.text = element_text(size =10),
#'         legend.title = element_text(size =10),
#'         panel.background = element_blank(),
#'         panel.grid.major = element_blank(), 
#'         panel.grid.minor = element_blank())+
#'   
#' 
#' ggsave(file="../../figs/scRNAseq/jpeg/fold_change_NT_INFg_all_cells_MHCI.jpeg", w=11, h=4)
#' ggsave(file="../../figs/scRNAseq/fold_change_NT_INFg_all_cells_MHCI.pdf", w=11, h=4)
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' #########################################################################################################
#' #Exploring cellular detection rate (CDR) from MAST
#' 
#' MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","MR1","CD1D","IRF1","NLRC5")
#' 
#' 
#' ## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
#' # but not normalized or scaled
#' d10x.stim<-readRDS("../../data/d10x_raw_merged.rds")
#' 
#' ## add cell type labels from split analysis
#' load(here("../../output/","celltype_labels_organoid.RData"))
#' identical(rownames(cell_typelabel),rownames(d10x.stim@meta.data))
#' cell_typelabel$cluster_ID[which(cell_typelabel$cluster_ID=="entero_stem")]<-"enterocyte_stem"
#' d10x.stim <- AddMetaData(d10x.stim, metadata = cell_typelabel)
#' 
#' ##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
#' # This is log(TP10K+1)
#' d10x.stim <- NormalizeData(d10x.stim,scale.factor = 10000, normalization.method = "LogNormalize")
#' 
#' 
#' ## testing factor
#' d10x.stim$cell_stim<-paste(d10x.stim$cluster_ID, d10x.stim$orig.ident, sep = "_")
#' Idents(d10x.stim) <- "cell_stim"
#' 
#' table(d10x.stim$cluster_ID, d10x.stim$orig.ident)
#' 
#' 
#' 
#' #MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene, 
#' #consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#' #and linear regression for the continuous process (i.e., the expression level). T
#' 
#' contrasts_celltype_stim<-rbind(
#'   combinations(n = 3, r = 2, v = unique(d10x.stim$cell_stim)[grep("crypt",unique(d10x.stim$cell_stim))], repeats.allowed = FALSE),
#'   combinations(n = 3, r = 2, v = unique(d10x.stim$cell_stim)[grep("TA",unique(d10x.stim$cell_stim))], repeats.allowed = FALSE),
#'   combinations(n = 3, r = 2, v = unique(d10x.stim$cell_stim)[grep("enterocyte_T|enterocyte_N|enterocyte_I",unique(d10x.stim$cell_stim))], repeats.allowed = FALSE),
#'   combinations(n = 3, r = 2, v = unique(d10x.stim$cell_stim)[grep("enterocyte_mt",unique(d10x.stim$cell_stim))], repeats.allowed = FALSE))
#' 
#' 
#' # this is 258,588 tests across all comparisons (21549 genes, 12 comparisons)
#' 
#' 
#' #############
#' ## Diff without CDR latent variable (nFeature)
#' #############
#' diff_exp_all<-lapply(1:nrow(contrasts_celltype_stim), function(x){
#'   de<-FindMarkers(d10x.stim, ident.1 = contrasts_celltype_stim[x,1], ident.2 = contrasts_celltype_stim[x,2], test.use = "MAST", verbose=F)
#'   print(paste(contrasts_celltype_stim[x,1],"vs", contrasts_celltype_stim[x,2],":", nrow(de), sep=" "))
#'   de$gene<-rownames(de)
#'   rownames(de)<-NULL
#'   de<-de[,c(6,1:5)]
#'   de$cell.1<-contrasts_celltype_stim[x,1]
#'   de$cell.2<-contrasts_celltype_stim[x,2]
#'   de})
#' 
#' 
#' diff_exp_all<-do.call(rbind, diff_exp_all)
#' 
#' 
#' table(diff_exp_all$gene)[which(table(diff_exp_all$gene)>10)]
#' 
#' 
#' #############
#' ## Diff with CDR latent variable (nFeature)
#' #############
#' diff_exp_all_latentvar<-lapply(1:nrow(contrasts_celltype_stim), function(x){
#'   de<-FindMarkers(d10x.stim, ident.1 = contrasts_celltype_stim[x,1], ident.2 = contrasts_celltype_stim[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
#'   print(paste(contrasts_celltype_stim[x,1],"vs", contrasts_celltype_stim[x,2],":", nrow(de), sep=" "))
#'   de$gene<-rownames(de)
#'   rownames(de)<-NULL
#'   de<-de[,c(6,1:5)]
#'   de$cell.1<-contrasts_celltype_stim[x,1]
#'   de$cell.2<-contrasts_celltype_stim[x,2]
#'   de})
#' 
#' 
#' diff_exp_all_latentvar<-do.call(rbind, diff_exp_all_latentvar)
#' 
#' 
#' table(diff_exp_all_latentvar$gene)[which(table(diff_exp_all_latentvar$gene)>10)]
#' 
#' diff_exp_all[grep("B2M",diff_exp_all$gene),]
#' diff_exp_all_latentvar[grep("B2M",diff_exp_all_latentvar$gene),]
#' 
#' diff_exp_all[grep("NLRC5",diff_exp_all$gene),]
#' diff_exp_all_latentvar[grep("NLRC5",diff_exp_all_latentvar$gene),]
#' 
#' 
#' check1<-diff_exp_all[which(diff_exp_all$cell.1=="crypt_IFNg" & diff_exp_all$cell.2=="crypt_TNFa"),]
#' check2<-diff_exp_all_latentvar[which(diff_exp_all_latentvar$cell.1=="crypt_IFNg" & diff_exp_all_latentvar$cell.2=="crypt_TNFa"),]
#' length(intersect(check2$gene, check1$gene))
#' 
#' check2<-check2[match(check1$gene, check2$gene),]
#' identical(check2$gene, check1$gene)
#' cor(check2$p_val,check2$p_val)
#' 
#' plt<-data.frame(gene=check1$gene,p_nolatent=check1$p_val, p_latent=check2$p_val)
#' ggplot(plt, aes(p_nolatent,p_latent ))+geom_point()
#' 
#' 
#' 
#' 
#' ### venn
#' crypt_NT_IFNg<-diff_exp_all_latentvar[which(diff_exp_all_latentvar$cell.1=="crypt_IFNg" & diff_exp_all_latentvar$cell.2=="crypt_NT"),]
#' TA_NT_IFNg<-diff_exp_all_latentvar[which(diff_exp_all_latentvar$cell.1=="TA_IFNg" & diff_exp_all_latentvar$cell.2=="TA_NT"),]
#' entero_NT_IFNg<-diff_exp_all_latentvar[which(diff_exp_all_latentvar$cell.1=="enterocyte_IFNg" & diff_exp_all_latentvar$cell.2=="enterocyte_NT"),]
#' enteromt_NT_IFNg<-diff_exp_all_latentvar[which(diff_exp_all_latentvar$cell.1=="enterocyte_mt_IFNg" & diff_exp_all_latentvar$cell.2=="enterocyte_mt_NT"),]
#' 
#' 
#' intersect(crypt_NT_IFNg$gene, TA_NT_IFNg$gene)
#' intersect(crypt_NT_IFNg$gene, entero_NT_IFNg$gene)
#' intersect(crypt_NT_IFNg$gene, enteromt_NT_IFNg$gene)
#' 
#' intersect(enteromt_NT_IFNg$gene, entero_NT_IFNg$gene)
#' 
#' diff_inallcell<-intersect(intersect(TA_NT_IFNg$gene, entero_NT_IFNg$gene),intersect(crypt_NT_IFNg$gene, enteromt_NT_IFNg$gene))
#' intersect(MHCI,diff_inallcell)
#' 
#' 
#' crypt_specific<-crypt_NT_IFNg$gene[which(!(crypt_NT_IFNg$gene%in%entero_NT_IFNg$gene | crypt_NT_IFNg$gene%in%TA_NT_IFNg$gene | crypt_NT_IFNg$gene%in%enteromt_NT_IFNg$gene))]
#' entero_specific<-entero_NT_IFNg$gene[which(!(entero_NT_IFNg$gene%in%crypt_NT_IFNg$gene | entero_NT_IFNg$gene%in%TA_NT_IFNg$gene | entero_NT_IFNg$gene%in%enteromt_NT_IFNg$gene))]
#' enteromt_specific<-enteromt_NT_IFNg$gene[which(!(enteromt_NT_IFNg$gene%in%crypt_NT_IFNg$gene | enteromt_NT_IFNg$gene%in%TA_NT_IFNg$gene | enteromt_NT_IFNg$gene%in%entero_NT_IFNg$gene))]
#' TA_specific<-TA_NT_IFNg$gene[which(!(TA_NT_IFNg$gene%in%entero_NT_IFNg$gene | TA_NT_IFNg$gene%in%crypt_NT_IFNg$gene | TA_NT_IFNg$gene%in%enteromt_NT_IFNg$gene))]
#' 
#' entero_general<-entero_NT_IFNg$gene[which(!(entero_NT_IFNg$gene%in%crypt_NT_IFNg$gene | entero_NT_IFNg$gene%in%TA_NT_IFNg$gene))]
#' 
#' intersect(MHCI,crypt_specific)
#' intersect(MHCI,entero_specific)
#' intersect(MHCI,enteromt_specific)
#' intersect(MHCI,TA_specific)
#' 
#' intersect(MHCI,entero_general)
#' 
#' 
#' 
#' ###### whats about cell cycle, does it differ between cell types and conditions?
#' s.genes <- cc.genes$s.genes
#' g2m.genes <- cc.genes$g2m.genes
#' 
#' d10x.stim <- CellCycleScoring(d10x.stim, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
#' 
#' 
#' ggplot(d10x.stim@meta.data, aes(orig.ident, S.Score))+geom_boxplot()+facet_wrap(~cluster_ID)
#' ggplot(d10x.stim@meta.data, aes(orig.ident, G2M.Score))+geom_boxplot()+facet_wrap(~cluster_ID)
#' 
#' 
#' 
#' 
#' 
#' 
