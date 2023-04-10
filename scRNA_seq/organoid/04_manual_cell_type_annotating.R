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
library(here)

options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))



########################################
#' Unsupervised clustering
########################################
NT.cells.integrated<-readRDS(here("data","NT_normalized_integrated.rds"))
NT.cells.integrated


#' ## percent expressing in SCT assay
NT.cells.exp<-as.data.frame(NT.cells.integrated[["SCT"]]@data)
#common_genes<-sapply(1:nrow(NT.cells.exp[1:100,]), function(x) length(which(NT.cells.exp[x,]==0)))
common_genes<-names(which(apply(NT.cells.exp, 1, function(x) length(which(x==0)))/ncol(NT.cells.exp)<0.9))
length(common_genes)

#'## Differenital by cluster (default to SCT for DE between clusters)
DefaultAssay(object = NT.cells.integrated) <- "SCT"
NT.cells.integrated <- PrepSCTFindMarkers(NT.cells.integrated)
NT.cells.markers_granular <- FindAllMarkers(NT.cells.integrated, features=common_genes, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25, min.cells.group=50)
table(NT.cells.markers_granular$cluster)
save(NT.cells.markers_granular, file=here("output","NT_cells_markers_granular.Rdata"))


load(file=here("output","NT_cells_markers_granular.Rdata"))

#'## order by FC and take top 100
NT.cells.markers_granular.sig<-NT.cells.markers_granular[which(NT.cells.markers_granular$p_val_adj<1e-50),]
table(NT.cells.markers_granular.sig$cluster)

NT.cells.markers_granular.sig100<-do.call(rbind,lapply(unique(NT.cells.markers_granular$cluster), function(x){
  clust_marker<-NT.cells.markers_granular.sig[which(NT.cells.markers_granular.sig$cluster==x),]
  clust_marker<-clust_marker[order(abs(clust_marker$avg_log2FC)),]
  clust_marker[1:100,]
}))

NT.cells.markers_granular.sig10<-do.call(rbind,lapply(unique(NT.cells.markers_granular$cluster), function(x){
  clust_marker<-NT.cells.markers_granular.sig[which(NT.cells.markers_granular.sig$cluster==x),]
  clust_marker<-clust_marker[order(abs(clust_marker$avg_log2FC)),]
  clust_marker[1:10,]
}))



rect<-ggplot()+geom_rect(aes(xmin=1:nrow(NT.cells.markers_granular.sig10), xmax=1:nrow(NT.cells.markers_granular.sig10)+1, ymin=1, ymax=2, 
                             fill=cluster),NT.cells.markers_granular.sig10, color="black") + geom_hline(yintercept=-1, color="white")+coord_flip()+
  theme_bw()+theme(legend.position = "none",
                   axis.title=element_blank(),
                   axis.text=element_blank(),
                   axis.ticks=element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_blank(),
                   panel.border = element_blank(),plot.margin = unit(c(-0.8,0,-0.4,0), "cm"))+
  geom_text(aes(x=seq(1,nrow(NT.cells.markers_granular.sig10),10)+5, y=0, label=unique(NT.cells.markers_granular.sig10$cluster)), size=4)



grid.arrange(
  rect,DotPlot(NT.cells.integrated, features = unique(NT.cells.markers_granular.sig10$gene))+coord_flip()+
    RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster"),
  widths=c(0.1,0.5))





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
                   panel.border = element_blank(),plot.margin = unit(c(0,5,0,0.75), "cm"))+
  geom_text(aes(x=seq(1,nrow(genes_df), 2)+1, y=0, label=genes_df$cellType[seq(1,nrow(genes_df), 2)]), size=4)+
  scale_fill_manual(values=c("#cf92beff","#364557","#9dbc26ff","#efe532ff","#F6A317","#89DFE2"))#


grid.arrange(
  DotPlot(NT.cells.integrated, features = genes_df$gene)+
    RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                            plot.margin = unit(c(0.3,0.5,0.1,0.1), "cm")),
  rect,
  heights=c(0.9,0.1))

ggsave(
  grid.arrange(
    DotPlot(NT.cells.integrated, features = genes_df$gene)+
      RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                              plot.margin = unit(c(0.3,0.5,0.1,0.1), "cm")),
    rect,
    heights=c(0.9,0.1)),file=here("figs","NT_markerdotplot_granular_SCT.pdf"), w=7,h=4)


ggsave(
  grid.arrange(
    DotPlot(NT.cells.integrated, features = genes_df$gene)+
      RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                              plot.margin = unit(c(0.3,0.5,0.1,0.1), "cm")),
    rect,
    heights=c(0.9,0.1)),file=here("figs/jpeg","NT_markerdotplot_granular_SCT.jpeg"), w=7,h=4)



#' Dev cell genes

genes_df_devCell<-data.frame(gene=c("EPCAM","FABP1","MKI67","LGR5","sct_ATOH1"),
                             cellType=c("All Epi","Enterocyte","Crypt Cell","Crypt Cell","Secretory Epi"))

NT.cells.integrated@assays$RNA@counts["ATOH1",1:10]

rect<-ggplot()+geom_rect(aes(xmin=1:nrow(genes_df_devCell), xmax=1:nrow(genes_df_devCell)+1, ymin=1, ymax=2, 
                             fill=cellType),genes_df_devCell, color="black") + geom_hline(yintercept=-1, color="white")+
  theme_bw()+theme(legend.position = "none",
                   axis.title=element_blank(),
                   axis.text=element_blank(),
                   axis.ticks=element_blank(),
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                   panel.background = element_blank(), axis.line = element_blank(),
                   panel.border = element_blank(),plot.margin = unit(c(0,5,0,0.75), "cm"))+
  geom_text(aes(x=seq(1,nrow(genes_df_devCell))+0.5, y=0, label=genes_df_devCell$cellType), size=4)+
  scale_fill_manual(values=c("cornflowerblue","#5e88be","#084594","#a6bddb"))#



grid.arrange(
  DotPlot(NT.cells.integrated, features = genes_df_devCell$gene)+
    RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                            plot.margin = unit(c(0.3,0.5,0.1,0.1), "cm")),
  rect,
  heights=c(0.9,0.1))

ggsave(
  grid.arrange(
    DotPlot(NT.cells.integrated, features = genes_df_devCell$gene)+
      RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                              plot.margin = unit(c(0.3,0.5,0.1,0.1), "cm")),
    rect,
    heights=c(0.9,0.1)), file=here("figs","NT_DevCellmarkerdotplot_granular_SCT.pdf"), w=7,h=4)

ggsave(
  grid.arrange(
    DotPlot(NT.cells.integrated, features = genes_df_devCell$gene)+
      RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                              plot.margin = unit(c(0.3,0.5,0.1,0.1), "cm")),
    rect,
    heights=c(0.9,0.1)), file=here("figs/jpeg","NT_DevCellmarkerdotplot_granular_SCT.jpeg"), w=7,h=4)



#'## smilie markers
#' 
#top markers
NT.cells.integrated@assays$RNA@counts["ATOH1",1:10]

## compare to markers from https://pubmed.ncbi.nlm.nih.gov/31348891/
epi.marker<-read.csv(here("data","NIHMS1532849-supplement-9.csv"),  header=T) 
#epi.marker_measured<-epi.marker[which(epi.marker$gene%in%rownames(NT.cells.integrated@assays$RNA@counts)),]
epi.marker_measured<-epi.marker[which(epi.marker$gene%in%NT.cells.markers_granular.sig100$gene),]



epi.marker_top<-do.call(rbind,lapply(1:length(unique(epi.marker_measured$ident)), function(y){
  celltype<-unique(epi.marker_measured$ident)[y]
  cell_genes<-head(epi.marker_measured[which(epi.marker_measured$ident==celltype),c("ident","gene")])
}))

epi.marker_top_plt<-epi.marker_top
epi.marker_top_plt$gene_order <- factor(epi.marker_top_plt$gene, levels = unique(epi.marker_top_plt$gene))


rect<-ggplot(epi.marker_top_plt, aes(ident, gene_order)) +
  geom_tile(aes(fill = ident), colour = "black")+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "none",
        axis.title=element_blank(),
        axis.ticks=element_blank(),
        panel.background = element_blank(), 
        panel.border = element_blank(),plot.margin = unit(c(1,1,1,1), "cm"))



grid.arrange(
  rect,
  DotPlot(NT.cells.integrated, features = unique(epi.marker_top$gene))+coord_flip()+
    RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(panel.border = element_blank(),plot.margin = unit(c(1,1,4,1), "cm")),
  widths=c(0.5,0.5))


ggsave(
  grid.arrange(
    rect,
    DotPlot(NT.cells.integrated, features = unique(epi.marker_top$gene))+coord_flip()+
      RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(panel.border = element_blank(),plot.margin = unit(c(1,1,4,1), "cm")),
    widths=c(0.3,0.4)), file=here("figs","NT_Smiliemarkerdotplot_granular_SCT.pdf"), w=15,h=15)

ggsave(
  grid.arrange(
    rect,
    DotPlot(NT.cells.integrated, features = unique(epi.marker_top$gene))+coord_flip()+
      RotatedAxis()+theme(axis.title.x=element_blank())+ylab("Cluster")+theme(panel.border = element_blank(),plot.margin = unit(c(1,1,4,1), "cm")),
    widths=c(0.3,0.4)), file=here("figs/jpeg","NT_Smiliemarkerdotplot_granular_SCT.jpeg"), w=15,h=15)







## significance of markers
NT.cells.markers_granular.sig[which(NT.cells.markers_granular.sig$gene%in%c(genes_df$gene, genes_df_devCell$gene)),]


#######################
#You can use the corrected log-normalized counts for differential expression and integration.
#However, in principle, it would be most optimal to perform these calculations directly on the
#residuals (stored in the scale.data slot) themselves. This is not currently supported in Seurat v3, but will be soon.

#log-normalized versions of corrected counts
NT.cells.exp<-as.data.frame(NT.cells.integrated[["SCT"]]@data)
table(unlist(NT.cells.exp["LGR5",]))

genes<-c(unique(NT.cells.markers_granular.sig10$gene),genes_df$gene, genes_df_devCell$gene)

NT.cells.exp.GOI<-NT.cells.exp[genes,]
NT.cells.exp.GOI$gene<-rownames(NT.cells.exp.GOI)
NT.cells.exp.GOI<-melt(NT.cells.exp.GOI)#

umap_mat<-as.data.frame(Embeddings(object = NT.cells.integrated, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-NT.cells.integrated@meta.data
meta$cell<-rownames(meta)

plt<-merge(NT.cells.exp.GOI, meta,by.x="variable", by.y="cell")
plt<-merge(plt, umap_mat,by.x="variable", by.y="cell")

# plt$color<-plt$value
# plt$color[which(plt$color==0)]<-NA

gene_plot<-function(gene_name){
  plt_gene<-plt[which(plt$gene==gene_name),]
  plt_gene<-plt_gene[order(plt_gene$value),]
  ggplot(plt_gene, aes(UMAP_1,UMAP_2, color=value))+
    geom_point(size=1.5)+theme_bw()+
    scale_colour_gradient2( low = "#2b8cbe",#blue
                            mid = "#f0f0f0", # grey
                            high = "#eb8423", #oragne
                            midpoint = max(plt_gene[,"value"])/2,
                            na.value = "#226d94",
                            name="Count")+theme_classic()+
    ggtitle(gene_name)+th+theme(plot.title = element_text(size = 16, face = "bold"))
}


tsne_mat<-as.data.frame(Embeddings(object = NT.cells.integrated, reduction = "tsne"))#
tsne_mat$cell<-rownames(tsne_mat)

plt_tsne<-merge(NT.cells.exp.GOI, meta,by.x="variable", by.y="cell")
plt_tsne<-merge(plt_tsne, tsne_mat,by.x="variable", by.y="cell")

gene_plot_tsne<-function(gene_name){
  plt_gene<-plt_tsne[which(plt_tsne$gene==gene_name),]
  plt_gene<-plt_gene[order(plt_gene$value),]
  ggplot(plt_gene, aes(tSNE_1,tSNE_2, color=value))+
    geom_point(size=1.5)+theme_bw()+
    scale_colour_gradient2( low = "#2b8cbe",#blue
                            mid = "#f0f0f0", # grey
                            high = "#eb8423", #oragne
                            midpoint = max(plt_gene[,"value"])/2,
                            na.value = "#226d94",
                            name="Count")+theme_classic()+
    ggtitle(gene_name)+th+theme(plot.title = element_text(size = 16, face = "bold"))
}


grid.arrange(gene_plot("LGR5"),gene_plot("ASCL2"),gene_plot("FABP1"),gene_plot("KRT19"),gene_plot("PCNA"),gene_plot("MKI67"),gene_plot("HELLS"),gene_plot("TOP2A"),ncol=2)

ggsave(file=here("figs","NT_markers.pdf"),
       grid.arrange(gene_plot("LGR5"),gene_plot("ASCL2"),gene_plot("FABP1"),gene_plot("KRT19"),gene_plot("PCNA"),gene_plot("MKI67"),gene_plot("HELLS"),gene_plot("TOP2A"),ncol=2),
       width=8,height=24)
ggsave(file=here("figs/jpeg","NT_markers.jpeg"),
       grid.arrange(gene_plot("LGR5"),gene_plot("ASCL2"),gene_plot("FABP1"),gene_plot("KRT19"),gene_plot("PCNA"),gene_plot("MKI67"),gene_plot("HELLS"),gene_plot("TOP2A"),ncol=2),
       width=8,height=16)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

cols = muted(gg_color_hue(11),l = 90, c = 30)


lapply(0:5, function(clust){
  markers<-rbind(NT.cells.markers_granular.sig10[which(NT.cells.markers_granular.sig10$cluster==clust),],NT.cells.markers_granular.sig[which(NT.cells.markers_granular.sig$gene%in%c(genes_df$gene, genes_df_devCell$gene)&NT.cells.markers_granular.sig$cluster==clust),])
  #write.table(markers$gene, file=paste(here("../../output/"),"top_markers_merged_cluster_",clust,".txt", sep=""),row.names = F, col.names = F, quote = F)
  
  write.table(paste(markers$gene, collapse=", "), file=paste(here("output/"),"top_markers_granular_cluster_",clust,".txt", sep=""),row.names = F, col.names = F, quote = F)
  
  key_marker<-NT.cells.markers_granular.sig[which(NT.cells.markers_granular.sig$gene%in%c(genes_df$gene, genes_df_devCell$gene)&NT.cells.markers_granular.sig$cluster==clust),]$gene
  top_cluster<-NT.cells.markers_granular.sig10[which(NT.cells.markers_granular.sig10$cluster==clust),]
  top_cluster[order(top_cluster$pct.2),]
  
  
  if(length(key_marker)==4){genes_to_plot<-key_marker}else{genes_to_plot<-c(key_marker,top_cluster[order(top_cluster$pct.2),"gene"][1:(4-length(key_marker))])}
  
  ggsave(do.call(grid.arrange,c(lapply(genes_to_plot,gene_plot),ncol=2)),
         file=paste(here("figs"),"/Cluster_",clust,"_granular_markers.pdf", sep=""), w=12,h=10)
  
  
  plt_clust<-plt[,c("variable","seurat_clusters","UMAP_1","UMAP_2")]
  plt_clust<-plt_clust[!duplicated(plt_clust),]
  
  cols[clust+1] = gg_color_hue(11)[clust+1]
  
  ggplot(plt_clust, aes(UMAP_1,UMAP_2, color=seurat_clusters))+
    geom_point(size=1.5, shape=20)+theme_bw()+
    theme_classic()+
    scale_color_manual(values=cols, name="Cluster")+th+theme(plot.title = element_text(size = 16, face = "bold"))
  
  ggsave(file=paste(here("figs"),"/Cluster_",clust,"_UMAP_granular.pdf", sep=""), w=6,h=5)})












##################################
#'# Cytocine treatments
##################################
genes_df<-data.frame(gene=c("LGR5","ASCL2","FABP1","KRT19","HELLS","PCNA","MKI67","TOP2A", "FCGBP","SPINK4", "SOX4","CD24"),
                     cellType=c("Stem","Stem","Early E","Early E","TA1","TA1","TA2","TA2","Goblet","Goblet","Tuft","Tuft"))

genes_df_devCell<-data.frame(gene=c("EPCAM","FABP1","MKI67","LGR5","sct_ATOH1"),
                             cellType=c("All Epi","Enterocyte","Crypt Cell","Crypt Cell","Secretory Epi"))



IFNg.cells.integrated<-readRDS(here("data","IFNg_normalized_integrated.rds"))

IFNg.cells.exp<-as.data.frame(IFNg.cells.integrated[["SCT"]]@data)
# ## percent expressing in SCT assay
common_genes<-names(which(apply(IFNg.cells.exp, 1, function(x) length(which(x==0)))/ncol(IFNg.cells.exp)<0.9))
length(common_genes)

## Differenital by cluster (default to SCT for DE between clusters)
DefaultAssay(object = IFNg.cells.integrated) <- "SCT"
IFNg.cells.integrated <- PrepSCTFindMarkers(IFNg.cells.integrated)
IFNg.cells.markers <- FindAllMarkers(IFNg.cells.integrated, features=common_genes, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25, min.cells.group=50)
table(IFNg.cells.markers$cluster)

# order by FC and take top 100
IFNg.cells.markers.sig<-IFNg.cells.markers[which(IFNg.cells.markers$p_val_adj<1e-50),]
table(IFNg.cells.markers.sig$cluster)

IFNg.cells.markers.sig100<-do.call(rbind,lapply(0:5, function(x){
  clust_marker<-IFNg.cells.markers.sig[which(IFNg.cells.markers.sig$cluster==x),]
  clust_marker<-clust_marker[order(abs(clust_marker$avg_log2FC)),]
  clust_marker[1:100,]
}))

IFNg.cells.markers.sig10<-do.call(rbind,lapply(0:5, function(x){
  clust_marker<-IFNg.cells.markers.sig[which(IFNg.cells.markers.sig$cluster==x),]
  clust_marker<-clust_marker[order(abs(clust_marker$avg_log2FC)),]
  clust_marker[1:10,]
}))


genes<-c(unique(IFNg.cells.markers.sig10$gene),genes_df$gene, genes_df_devCell$gene)

IFNg.cells.exp.GOI<-IFNg.cells.exp[genes,]
IFNg.cells.exp.GOI$gene<-rownames(IFNg.cells.exp.GOI)
IFNg.cells.exp.GOI<-melt(IFNg.cells.exp.GOI)#

umap_mat<-as.data.frame(Embeddings(object = IFNg.cells.integrated, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-IFNg.cells.integrated@meta.data
meta$cell<-rownames(meta)

plt<-merge(IFNg.cells.exp.GOI, meta,by.x="variable", by.y="cell")
plt<-merge(plt, umap_mat,by.x="variable", by.y="cell")

gene_plot<-function(gene_name){
  plt_gene<-plt[which(plt$gene==gene_name),]
  plt_gene<-plt_gene[order(plt_gene$value),]
  ggplot(plt_gene, aes(UMAP_1,UMAP_2, color=value))+
    geom_point(size=1.5)+theme_bw()+
    scale_colour_gradient2( low = "#2b8cbe",#blue
                            mid = "#f0f0f0", # grey
                            high = "#eb8423", #oragne
                            midpoint = max(plt_gene[,"value"])/2,
                            na.value = "#226d94",
                            name="Count")+theme_classic()+
    ggtitle(gene_name)+th+theme(plot.title = element_text(size = 16, face = "bold"))
}

grid.arrange(gene_plot("LGR5"),gene_plot("ASCL2"),gene_plot("FABP1"),gene_plot("KRT19"),gene_plot("PCNA"),gene_plot("MKI67"),gene_plot("HELLS"),gene_plot("TOP2A"),ncol=2)

ggsave(file=here("figs","IFNg_markers.pdf"),
       grid.arrange(gene_plot("LGR5"),gene_plot("ASCL2"),gene_plot("FABP1"),gene_plot("KRT19"),gene_plot("PCNA"),gene_plot("MKI67"),gene_plot("HELLS"),gene_plot("TOP2A"),ncol=2),
       width=8,height=24)
ggsave(file=here("figs/jpeg","IFNg_markers.jpeg"),
       grid.arrange(gene_plot("LGR5"),gene_plot("ASCL2"),gene_plot("FABP1"),gene_plot("KRT19"),gene_plot("PCNA"),gene_plot("MKI67"),gene_plot("HELLS"),gene_plot("TOP2A"),ncol=2),
       width=8,height=16)



#############

TNFa.cells.integrated<-readRDS(here("data","TNFa_normalized_integrated.rds"))

TNFa.cells.exp<-as.data.frame(TNFa.cells.integrated[["SCT"]]@data)
# ## percent expressing in SCT assay
common_genes<-names(which(apply(TNFa.cells.exp, 1, function(x) length(which(x==0)))/ncol(TNFa.cells.exp)<0.9))
length(common_genes)

## Differenital by cluster (default to SCT for DE between clusters)
DefaultAssay(object = TNFa.cells.integrated) <- "SCT"
TNFa.cells.integrated <- PrepSCTFindMarkers(TNFa.cells.integrated)
TNFa.cells.markers <- FindAllMarkers(TNFa.cells.integrated, features=common_genes, only.pos=TRUE, min.pct=0.1, logfc.threshold=0.25, min.cells.group=50)
table(TNFa.cells.markers$cluster)

# order by FC and take top 100
TNFa.cells.markers.sig<-TNFa.cells.markers[which(TNFa.cells.markers$p_val_adj<1e-50),]
table(TNFa.cells.markers.sig$cluster)

TNFa.cells.markers.sig100<-do.call(rbind,lapply(0:5, function(x){
  clust_marker<-TNFa.cells.markers.sig[which(TNFa.cells.markers.sig$cluster==x),]
  clust_marker<-clust_marker[order(abs(clust_marker$avg_log2FC)),]
  clust_marker[1:100,]
}))

TNFa.cells.markers.sig10<-do.call(rbind,lapply(0:5, function(x){
  clust_marker<-TNFa.cells.markers.sig[which(TNFa.cells.markers.sig$cluster==x),]
  clust_marker<-clust_marker[order(abs(clust_marker$avg_log2FC)),]
  clust_marker[1:10,]
}))


genes<-c(unique(TNFa.cells.markers.sig10$gene),genes_df$gene, genes_df_devCell$gene)

TNFa.cells.exp.GOI<-TNFa.cells.exp[genes,]
TNFa.cells.exp.GOI$gene<-rownames(TNFa.cells.exp.GOI)
TNFa.cells.exp.GOI<-melt(TNFa.cells.exp.GOI)#

umap_mat<-as.data.frame(Embeddings(object = TNFa.cells.integrated, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-TNFa.cells.integrated@meta.data
meta$cell<-rownames(meta)

plt<-merge(TNFa.cells.exp.GOI, meta,by.x="variable", by.y="cell")
plt<-merge(plt, umap_mat,by.x="variable", by.y="cell")

gene_plot<-function(gene_name){
  plt_gene<-plt[which(plt$gene==gene_name),]
  plt_gene<-plt_gene[order(plt_gene$value),]
  ggplot(plt_gene, aes(UMAP_1,UMAP_2, color=value))+
    geom_point(size=1.5)+theme_bw()+
    scale_colour_gradient2( low = "#2b8cbe",#blue
                            mid = "#f0f0f0", # grey
                            high = "#eb8423", #oragne
                            midpoint = max(plt_gene[,"value"])/2,
                            na.value = "#226d94",
                            name="Count")+theme_classic()+
    ggtitle(gene_name)+th+theme(plot.title = element_text(size = 16, face = "bold"))
}

grid.arrange(gene_plot("LGR5"),gene_plot("ASCL2"),gene_plot("FABP1"),gene_plot("KRT19"),gene_plot("PCNA"),gene_plot("MKI67"),gene_plot("HELLS"),gene_plot("TOP2A"),ncol=2)

ggsave(file=here("figs","TNFa_markers.pdf"),
       grid.arrange(gene_plot("LGR5"),gene_plot("ASCL2"),gene_plot("FABP1"),gene_plot("KRT19"),gene_plot("PCNA"),gene_plot("MKI67"),gene_plot("HELLS"),gene_plot("TOP2A"),ncol=2),
       width=8,height=24)
ggsave(file=here("figs/jpeg","TNFa_markers.jpeg"),
       grid.arrange(gene_plot("LGR5"),gene_plot("ASCL2"),gene_plot("FABP1"),gene_plot("KRT19"),gene_plot("PCNA"),gene_plot("MKI67"),gene_plot("HELLS"),gene_plot("TOP2A"),ncol=2),
       width=8,height=16)

