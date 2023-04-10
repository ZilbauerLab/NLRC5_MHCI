library(dplyr)
library(Seurat)
library(patchwork)
library(here)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(colorspace)
library(gtools)

options(stringsAsFactors = FALSE)

#source(here("R_functions/","pretty_plots.R"))
source(here("../EBI/scRNAseq_codon/scripts/00_pretty_plots.R"))


adult_d10x.primary<-readRDS(here("data/adult_intestine","local.rds"))
adult_d10x.primary

print(head(adult_d10x.primary@meta.data))

print(table(adult_d10x.primary@meta.data$category))
print(table(adult_d10x.primary@meta.data$tissue))

adult_d10x_epi<-subset(adult_d10x.primary, subset = category %in% c("Epithelial"))
adult_d10x_epi<-subset(adult_d10x_epi, subset = tissue %in% c("large intestine","rectum","small intestine"))

print(table(adult_d10x_epi@meta.data$Age_group))

adult_d10x_epi<-subset(adult_d10x_epi, subset = Age_group %in% c("Adult"))
adult_d10x_epi

DimPlot(adult_d10x_epi, reduction = "umap", group.by = "tissue", pt.size=0.25, label=T)
ggsave(file=here("figs","adult_primary_UMAP_cluster.pdf"), w=6,h=5)
ggsave(file=here("figs/jpeg","adult_primary_UMAP_cluster.jpeg"), w=6,h=5)

head(colnames(adult_d10x_epi))
head(rownames(adult_d10x_epi))

MHCI = c('ENSG00000204642', 'ENSG00000204632', 'ENSG00000206503', 'ENSG00000204592', 
         'ENSG00000204525', 'ENSG00000234745',"ENSG00000168394","ENSG00000204267",
         "ENSG00000240065","ENSG00000204264","ENSG00000166710","ENSG00000125347","ENSG00000140853")

crypt_villis = c( "ENSG00000280501" , "ENSG00000111701", "ENSG00000136826", "ENSG00000204936", "ENSG00000154269",
                  "ENSG00000084674" ,"ENSG00000196839" ,"ENSG00000137860", "ENSG00000164120" ,"ENSG00000079385", "ENSG00000171431", "ENSG00000135318",
                  "ENSG00000165949" ,"ENSG00000259823" ,"ENSG00000145287" ,"ENSG00000015520" ,"ENSG00000135549" ,"ENSG00000110244", "ENSG00000118137",
                  "ENSG00000007306" ,"ENSG00000105388" ,"ENSG00000086548" ,"ENSG00000146648", "ENSG00000073737" ,"ENSG00000155366", "ENSG00000117472")


adult_d10x_epi <- AddModuleScore(
  object = adult_d10x_epi,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score'
)

adult_d10x_epi <- AddModuleScore(
  object = adult_d10x_epi,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)


adult_d10x_epi<-subset(adult_d10x_epi, subset = cell_type != "progenitor cell of endocrine pancreas")
adult_d10x_epi<-subset(adult_d10x_epi, subset = cell_type != "GIP cell")
adult_d10x_epi<-subset(adult_d10x_epi, subset = cell_type != "glial cell")

adult_d10x_epi <-FindVariableFeatures(adult_d10x_epi, selection.method = "vst", nfeatures = 2000)
adult_d10x_epi <-ScaleData(adult_d10x_epi)
adult_d10x_epi <- RunPCA(adult_d10x_epi, npcs = 30, verbose = FALSE)
adult_d10x_epi <- RunUMAP(adult_d10x_epi, reduction = "pca", dims = 1:30)

# DimPlot(adult_d10x_epi, reduction = "umap", group.by = "tissue", pt.size=0.25, label=T)
# ggsave(file=here("figs","adult_primary_epi_UMAP_tissue.pdf"), w=8,h=5)
# ggsave(file=here("figs/jpeg","adult_primary_UMAP_tissue.jpeg"), w=8,h=5)
# 
# DimPlot(adult_d10x_epi, reduction = "umap", group.by = "cell_type", pt.size=0.25, label=T)
# ggsave(file=here("figs","adult_primary_epi_UMAP_cell_type.pdf"), w=10,h=5)
# ggsave(file=here("figs/jpeg","adult_primary_UMAP_cell_type.jpeg"), w=10,h=5)
# 
# FeaturePlot(adult_d10x_epi, features = "MHCI_score1",reduction = "umap", min.cutoff = "q9", pt.size=1)
# ggsave(file=here("figs","adult_primary_epi_UMAP_MHCI.pdf"), w=6,h=5)
# ggsave(file=here("figs/jpeg","adult_primary_UMAP_MHCI.jpeg"), w=6,h=5)
# 
# 
# ggplot(adult_d10x_epi@meta.data, aes(reorder(cell_type, MHCI_score1, FUN = median),MHCI_score1))+
#   geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=tissue))+
#   theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),
#                       axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~tissue, scales="free_x")
# ggsave(file=here("figs/adult_primary_epi_MHCI_boxplot_celltype.pdf"), w=30, h=10)
# ggsave(file=here("figs/jpeg/adult_primary_epi_MHCI_boxplot_celltype.jpeg"), w=30, h=10)
# 
# 
# #### data frame for plotting
# umap_mat<-as.data.frame(Embeddings(object = adult_d10x_epi, reduction = "umap"))#
# umap_mat$cell<-rownames(umap_mat)
# 
# meta<-adult_d10x_epi@meta.data
# meta$cell<-rownames(meta)
# 
# plt_UMAP<-merge(meta, umap_mat, by="cell")
# 
# 
# ###### summary plot
# table(adult_d10x_epi$donor_id, adult_d10x_epi$tissue)
# adult_d10x_epi
# 
# plt<-adult_d10x_epi@meta.data[,c("MHCI_score1","crypt_villis_score1","author_cell_type","tissue")]
# save(plt,plt_UMAP, file=here("data/adult_intestine","plt_data.RData"))
# 
# load("../EBI/MHCI/more_single_cell/plt_data.RData")
# 
# 
# 
# 
# ##################################################################################################
# 
#             table(plt$cell_type)
#             table(as.character(plt_UMAP$author_cell_type), as.character(plt_UMAP$cell_type))
#             table(as.character(plt_UMAP$cell_type), as.character(plt_UMAP$author_cell_type))
# 
# plt_UMAP<-plt_UMAP[which(!(plt_UMAP$author_cell_type%in%c("M/X cells (MLN/GHRL+)","Progenitor (NEUROG3+)"))),]
# plt_UMAP$author_cell_type<-as.factor(as.character(plt_UMAP$author_cell_type))
# levels(plt_UMAP$author_cell_type)<-c("BEST2+ Goblet cell","BEST4+ epithelial","Colonocyte","Enteroendocrine","Enteroendocrine",
#                                      "Enteroendocrine","Enterocyte","Goblet cell","Enteroendocrine","Enteroendocrine", "Enteroendocrine", 
#                                      "Microfold cell","Enteroendocrine","Paneth","Stem cells","TA","Tuft")
# levels(plt_UMAP$author_cell_type)
# cols_manual<-c("#a30443ff", "#87d235ff","#04b026","#5c1e6cff","#208d43ff", "#ad0044ff","#5b1f6cff","#d9c025ff","#6493ebff","#ef6a13ff","#a03c9fff")
# names(cols_manual) <- levels(plt_UMAP$author_cell_type)
# fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
# colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)
# 
# ggplot(plt_UMAP, aes(UMAP_1,UMAP_2, color=author_cell_type))+
#   geom_point(size=1.5)+colcale_cols_manual+
#   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))
# ggsave(file=here("figs", "adult_primary_epi_celltype.pdf"), w=7,h=5)
# ggsave(file=here("figs/jpeg", "adult_primary_epi_celltype.jpeg"), w=7,h=5,bg="white")
# 
# ggplot(plt_UMAP, aes(UMAP_1,UMAP_2, color=tissue))+
#   geom_point(size=1.5)+scale_color_manual(values=c("#a6d96a","darkgreen","cornflowerblue"))+
#   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))
# ggsave(file=here("figs", "adult_primary_epi_tissue.pdf"), w=7,h=5)
# ggsave(file=here("figs/jpeg", "adult_primary_epi_tissue.jpeg"), w=7,h=5,bg="white")
# 
# ##################################################################################################
# 
# ## simplify cell type
# plt<-plt[which(!(plt$author_cell_type%in%c("M/X cells (MLN/GHRL+)","Progenitor (NEUROG3+)"))),]
# 
# ## merge all enteroendocrine types
# plt$author_cell_type<-as.factor(as.character(plt$author_cell_type))
# levels(plt$author_cell_type)<-c("BEST2+ Goblet cell","BEST4+ epithelial","Colonocyte","Enteroendocrine","Enteroendocrine",
#                                      "Enteroendocrine","Enterocyte","Goblet cell","Enteroendocrine","Enteroendocrine", "Enteroendocrine", 
#                                      "Microfold cell","Enteroendocrine","Paneth","Stem cells","TA","Tuft")
# plt$author_cell_type<-as.character(plt$author_cell_type)
# 
# plt$tissue<-as.factor(as.character(plt$tissue))
# levels(plt$tissue)
# levels(plt$tissue)<-c("Large\nIntestine","Rectum","Small\nIntestine")
# 
# scatter<-ggplot(plt, aes(crypt_villis_score1, MHCI_score1, fill=tissue, color=tissue))+geom_point(shape=21, color="black", size=1.5)+
#   theme_bw()+th+theme(legend.position = "none",plot.margin = unit(c(0,0,1.4,4.6), "cm"))+
#   scale_color_manual(values=c("#a6d96a","#034203","cornflowerblue"))+scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"))+  stat_smooth(method="lm",se=F)+
#   xlab("Crypt-Villus Score")+ylab("MHC I Score")
# 
# 
# violin_celltype<-ggplot(plt, aes(reorder(author_cell_type, crypt_villis_score1), crypt_villis_score1, fill=author_cell_type))+
#   geom_violin()+scale_x_discrete(position = "bottom") +coord_flip()+fillscale_cols_manual+
#   theme_bw()+th+ylab("")+xlab("")+theme(axis.title.x=element_blank(),
#                                         axis.text.x=element_blank(),
#                                         axis.ticks.x=element_blank(),
#                                         legend.position = "none",
#                                         plot.margin = unit(c(1,0,0,1.7), "cm"))
# 
# box_diagnosis<-ggplot(plt, aes(reorder(tissue,MHCI_score1, mean), MHCI_score1, fill=tissue))+geom_boxplot()+xlab("")+
#   scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), name="Segment")+theme_bw()+th+theme(axis.title.y =element_blank(),
#                                                                                                           axis.text.y=element_blank(),
#                                                                                                           axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#                                                                                                           axis.ticks.y=element_blank(),
#                                                                                                           plot.margin = unit(c(0,0.2,0.45,0), "cm"))
# 
# mn_cryptvillus<-tapply(plt$crypt_villis_score1, plt$author_cell_type, mean)
# cryptvillus_order<-names(mn_cryptvillus)[(order(mn_cryptvillus))]
# 
# plt_count<-melt(table(plt$tissue, plt$author_cell_type))
# colnames(plt_count)<-c("Segment","Cell_Type","count")
# plt_count$Cell_Type<-factor(plt_count$Cell_Type, cryptvillus_order)
# 
# placeholder<-ggplot(plt_count, aes(Segment, Cell_Type))+geom_tile(color="grey50", fill="white",size=0.4)+
#   geom_text(aes(label=count),size=2.75,color="grey40")+ggtitle(" Cell Number")+
#   theme_void()+theme(axis.title=element_blank(),
#                      axis.text=element_blank(),
#                      axis.ticks=element_blank(),
#                      plot.margin = unit(c(0.4,3.3,0,0), "cm"),
#                      plot.title = element_text(size = 10))
# placeholder
# 
# grid_fig<-plot_grid(violin_celltype, placeholder, scatter, box_diagnosis, 
#                     #align = 'hv',
#                     ncol = 2, axis="lr", 
#                     rel_heights = c(1,2),
#                     rel_widths = c(4,1.5))
#                             ##################################################################################################
# 
#                             ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_epi_MHCI_crpyt_summary.pdf"),grid_fig, w=10,h=10)
#                             ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_epi_MHCI_crpyt_summary.jpeg"),grid_fig, w=9,h=8,bg="white")
#                   
#                             ##################################################################################################
#                             
# ggsave(file=here("figs", "adult_primary_epi_MHCI_crpyt_summary.pdf"),grid_fig, w=9,h=8)
# ggsave(file=here("figs/jpeg", "adult_primary_epi_MHCI_crpyt_summary.jpeg"),grid_fig, w=9,h=8,bg="white")
#                   
# cor(plt$MHCI_score1, plt$crypt_villis_score1, method = "spearman")
# 
# 
# library(stats)
# summary(aov(plt$MHCI_score1 ~ plt$tissue))
# TukeyHSD(aov(plt$MHCI_score1 ~ plt$tissue))
# 
# 
# ggplot(plt, aes(tissue,MHCI_score1))+geom_violin(color="lightgrey",fill="lightgrey")+geom_boxplot(aes(fill=tissue))+facet_wrap(~author_cell_type)+scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), name="Segment")+theme_bw()+th+ylab("MHC I Score")+xlab("")
# ggsave(file=here("figs", "adult_primary_epi_MHCI_by_celltype.pdf"), w=9,h=8)
# ggsave(file=here("figs/jpeg", "adult_primary_epi_MHCI__celltype.jpeg"), w=9,h=8,bg="white")
# 




##########
## Scatter spit by segment
##########
load("../EBI/MHCI/more_single_cell/plt_data.RData")

table(plt$cell_type)
table(as.character(plt_UMAP$author_cell_type), as.character(plt_UMAP$cell_type))
table(as.character(plt_UMAP$cell_type), as.character(plt_UMAP$author_cell_type))

plt_UMAP<-plt_UMAP[which(!(plt_UMAP$author_cell_type%in%c("M/X cells (MLN/GHRL+)","Progenitor (NEUROG3+)"))),]
plt_UMAP$author_cell_type<-as.factor(as.character(plt_UMAP$author_cell_type))
levels(plt_UMAP$author_cell_type)<-c("BEST2+ Goblet cell","BEST4+ epithelial","Colonocyte","Enteroendocrine","Enteroendocrine",
                                     "Enteroendocrine","Enterocyte","Goblet cell","Enteroendocrine","Enteroendocrine", "Enteroendocrine",
                                     "Microfold cell","Enteroendocrine","Paneth","Stem cells","TA","Tuft")
levels(plt_UMAP$author_cell_type)
cols_manual<-c("#a30443ff", "#87d235ff","#04b026","#5c1e6cff","#208d43ff", "#ad0044ff","#5b1f6cff","#d9c025ff","#6493ebff","#ef6a13ff","#a03c9fff")
names(cols_manual) <- levels(plt_UMAP$author_cell_type)
fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)



## simplify cell type
plt<-plt[which(!(plt$author_cell_type%in%c("M/X cells (MLN/GHRL+)","Progenitor (NEUROG3+)"))),]

## merge all enteroendocrine types
plt$author_cell_type<-as.factor(as.character(plt$author_cell_type))
levels(plt$author_cell_type)<-c("BEST2+ Goblet cell","BEST4+ epithelial","Colonocyte","Enteroendocrine","Enteroendocrine",
                                "Enteroendocrine","Enterocyte","Goblet cell","Enteroendocrine","Enteroendocrine", "Enteroendocrine",
                                "Microfold cell","Enteroendocrine","Paneth","Stem cells","TA","Tuft")
plt$author_cell_type<-as.character(plt$author_cell_type)

plt$tissue<-as.factor(as.character(plt$tissue))
levels(plt$tissue)
levels(plt$tissue)<-c("Large\nIntestine","Rectum","Small\nIntestine")

plt$author_cell_type<-factor(plt$author_cell_type, levels=rev(c("BEST4+ epithelial","Colonocyte","Enterocyte","BEST2+ Goblet cell","Tuft","Enteroendocrine",
                                                            "Microfold cell","Goblet cell","Paneth","TA","Stem cells")))




seg_plt<-function(segment, colour){
  plt_sgement<-plt[which(plt$tissue==segment),]
  
  scatter<-ggplot(plt_sgement, aes(crypt_villis_score1, MHCI_score1, fill=tissue, color=tissue))+geom_point(shape=21, color="black", size=1.5)+
    theme_bw()+th+theme(legend.position = "none",plot.margin = unit(c(0,1,1,0), "cm"))+
    scale_color_manual(values="grey70")+scale_fill_manual(values=colour)+  stat_smooth(method="lm",se=F)+
    xlab("Crypt-Villus Score")+ylab("MHC I Score")+
    xlim(c(round(min(plt$crypt_villis_score1)), round(max(plt$crypt_villis_score1))))+
    ylim(c(round(min(plt$MHCI_score1)), round(max(plt$MHCI_score1))))+
    annotate("text", x=1.1, y=-0.75, label=paste("Rs = ",signif(cor.test(plt_sgement$crypt_villis_score1, plt_sgement$MHCI_score1, method="spearman")$estimate,2),";",
                                                 "p = ",signif(cor.test(plt_sgement$crypt_villis_score1, plt_sgement$MHCI_score1, method="spearman")$p.value,1)))
  
  plt_count<-melt(table( plt_sgement$author_cell_type))
  colnames(plt_count)<-c("author_cell_type","count")
  plt_count$type_count<-paste(plt_count$author_cell_type, "\n",comma(plt_count$count), sep="")
  
  plt_sgement$author_cell_type_count<- plt_sgement$author_cell_type
  levels(plt_sgement$author_cell_type_count)<- plt_count$type_count
  
  violin_celltype<-ggplot(plt_sgement, aes(author_cell_type_count, crypt_villis_score1, fill=author_cell_type))+
    geom_violin()+scale_x_discrete(position = "bottom") +coord_flip()+fillscale_cols_manual+
    theme_bw()+th+ylab("")+xlab("")+theme(axis.title.x=element_blank(),
                                          axis.text.x=element_blank(),
                                          axis.ticks.x=element_blank(),
                                          legend.position = "none",
                                          plot.margin = unit(c(1,1,1,0), "cm"))#+    geom_text(data=plt_count, aes(y=max(plt_sgement$crypt_villis_score1)*0.8, x=author_cell_type, label=comma(count)), vjust=-0.5)
  violin_celltype
  
  grid_fig<-plot_grid(violin_celltype, scatter,
                      align = 'v',
                      ncol = 1, axis="lr",
                      rel_heights = c(1,1),
                      rel_widths = c(4,1.5))
  grid_fig}

seg_plt("Large\nIntestine", "#a6d96a")

split_scatter<-plot_grid(
          seg_plt("Rectum", "darkgreen"),
          seg_plt("Large\nIntestine", "#a6d96a"),
          seg_plt("Small\nIntestine","cornflowerblue"), ncol=3)

split_scatter


ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_epi_MHCI_crpyt_summary_split.pdf"),split_scatter, w=15,h=8)
ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_epi_MHCI_crpyt_summary_split.jpeg"),split_scatter, w=15,h=8,bg="white")



comp_simple<-list(c("Large\nIntestine","Rectum"),c("Large\nIntestine","Small\nIntestine"),c("Small\nIntestine","Rectum"))
box_segment<-ggplot(plt, aes(reorder(tissue,MHCI_score1, mean), MHCI_score1, fill=tissue))+geom_boxplot()+xlab("")+ylab("MHC I Score")+
  scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), name="Segment")+theme_bw()+th+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
                                                                                                          legend.position = "none")+
  geom_signif(comparisons = comp_simple, step_increase = 0.05,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")

ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_epi_MHCI_crpyt_summary_box_segment.pdf"),box_segment, w=2,h=5)
ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_epi_MHCI_crpyt_summary_box_segment.jpeg"),box_segment, w=2,h=5,bg="white")



library(stats)
summary(aov(plt$MHCI_score1 ~ plt$tissue))
TukeyHSD(aov(plt$MHCI_score1 ~ plt$tissue))


ggplot(plt, aes(tissue,MHCI_score1))+geom_violin(color="lightgrey",fill="lightgrey")+geom_boxplot(aes(fill=tissue))+facet_wrap(~author_cell_type)+scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), name="Segment")+theme_bw()+th+ylab("MHC I Score")+xlab("")
ggsave(file=here("figs", "adult_primary_epi_MHCI_by_celltype.pdf"), w=9,h=8)
ggsave(file=here("figs/jpeg", "adult_primary_epi_MHCI__celltype.jpeg"), w=9,h=8,bg="white")



##### how many samples
load(here("/home/redgar/Documents/EBI/MHCI/more_single_cell","plt_MHCIexp_data.RData"))
meta<-plt_exp[,c("donor_id","tissue")]
meta<-meta[!duplicated(meta),]

table(meta$tissue)
length(unique(meta$donor_id))

###############
## Differential expression and UMAP epithelial only individual genes
###############
d10x.exp.GOI<-FetchData(object = adult_d10x_epi, vars = MHCI)
d10x.exp.GOI$gene<-rownames(d10x.exp.GOI)
d10x.exp.GOI<-melt(d10x.exp.GOI)#
head(d10x.exp.GOI)

umap_mat<-as.data.frame(Embeddings(object = adult_d10x_epi, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-adult_d10x_epi@meta.data
meta$cell<-rownames(meta)

plt<-merge(meta, umap_mat, by="cell")

plt_exp<-merge(d10x.exp.GOI, plt,by.x="variable", by.y="cell")
save(plt_exp, file=here("data/adult_intestine","plt_MHCIexp_data.RData"))


load(here("/home/redgar/Documents/EBI/MHCI/more_single_cell","plt_MHCIexp_data.RData"))

plt_exp$tissue<-as.factor(as.character(plt_exp$tissue))
levels(plt_exp$tissue)
levels(plt_exp$tissue)<-c("Large\nIntestine","Rectum","Small\nIntestine")

MHCI_df<-data.frame(ens_gene=MHCI, gene_name=c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5"))
plt_exp<-merge(plt_exp, MHCI_df, by.x="gene", by.y="ens_gene")


plt_gene_UMAP<-function(gene){
  plt_gene<-plt_exp[which(plt_exp$gene_name==gene),]

  boxplot<-ggplot(plt_gene, aes(tissue, value))+geom_violin(fill="grey70", color="grey70")+geom_boxplot(aes(fill=tissue),width=0.2,outlier.shape = NA)+th_present+theme_bw()+
    scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), guide="none")+ylab("Count")+xlab("Segment")
  
  plt_gene$value[which(plt_gene$value==0)]<-NA
  plt_gene<-plt_gene[rev(order(plt_gene$value)),]
  UMAP_gene<-ggplot(plt_gene, aes(UMAP_1, UMAP_2, color=value))+geom_point()+
    th_present+theme_bw()+scale_color_continuous_sequential("ag_GrnYl",begin=0.0001, rev = F,na.value = "grey80")+
    xlab("UMAP 1")+ylab("UMAP 2")
  
  UMAP_segemnt<-ggplot(plt_gene, aes(UMAP_1, UMAP_2, color=tissue))+geom_point()+
    th_present+theme_bw()+scale_color_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), guide="none")+
    xlab("UMAP 1")+ylab("UMAP 2")
  
  plot_grid(UMAP_segemnt,UMAP_gene,boxplot, ncol=3, rel_widths = c(0.8,1,0.5))}

plt_gene_UMAP("NLRC5")
ggsave(file=here("/home/redgar/Documents/EBI/MHCI/more_single_cell", "adult_primary_epi_MHCI_celltype_NLRC5.jpeg"), w=9,h=3,bg="white")

plt_gene_UMAP("B2M")
ggsave(file=here("/home/redgar/Documents/EBI/MHCI/more_single_cell", "adult_primary_epi_MHCI_celltype_B2M.jpeg"), w=9,h=3,bg="white")

#########
## DE
#########
## simplify cell type
d10x@meta.data<-d10x@meta.data[which(!(d10x@meta.data$author_cell_type%in%c("M/X cells (MLN/GHRL+)","Progenitor (NEUROG3+)"))),]

## merge all enteroendocrine types
d10x@meta.data$author_cell_type<-as.factor(as.character(d10x@meta.data$author_cell_type))
levels(d10x@meta.data$author_cell_type)<-c("BEST2+ Goblet cell","BEST4+ epithelial","Colonocyte","Enteroendocrine","Enteroendocrine",
                                     "Enteroendocrine","Enterocyte","Goblet cell","Enteroendocrine","Enteroendocrine", "Enteroendocrine",
                                     "Microfold cell","Enteroendocrine","Paneth","Stem cells","TA","Tuft")
d10x@meta.data$author_cell_type<-as.character(d10x@meta.data$author_cell_type)

d10x@meta.data$tissue<-as.factor(as.character(d10x@meta.data$tissue))
levels(d10x@meta.data$tissue)
levels(d10x@meta.data$tissue)<-c("Large\nIntestine","Rectum","Small\nIntestine")

## testing factor
Idents(adult_d10x_epi) <- "tissue"

table(adult_d10x_epi$author_cell_type, adult_d10x_epi$tissue)


#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 

cell_types<-unique(adult_d10x_epi$author_cell_type)

contrasts_celltype_tissue<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = adult_d10x_epi$tissue[grep(cell_types[x],adult_d10x_epi$tissue)], repeats.allowed = FALSE)}))

contrasts_celltype_tissue

nrow(contrasts_celltype_tissue)

diff_exp_all<-lapply(1:nrow(contrasts_celltype_tissue), function(x){
  de<-FindMarkers(adult_d10x_epi, ident.1 = contrasts_celltype_tissue[x,1], ident.2 = contrasts_celltype_tissue[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  print(paste(contrasts_celltype_tissue[x,1],"vs", contrasts_celltype_tissue[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_tissue[x,1]
  de$cell.2<-contrasts_celltype_tissue[x,2]
  de})

head(diff_exp_all)