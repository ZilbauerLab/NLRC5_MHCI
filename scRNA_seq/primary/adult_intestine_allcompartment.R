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
# 
# print(head(adult_d10x.primary@meta.data))
# 
# print(table(adult_d10x.primary@meta.data$category))
# print(table(adult_d10x.primary@meta.data$tissue))
# print(table(adult_d10x.primary@meta.data$author_cell_type, adult_d10x.primary@meta.data$category))
# print(table(adult_d10x.primary@meta.data$author_cell_type, adult_d10x.primary@meta.data$tissue))
# 

adult_d10x.primary<-subset(adult_d10x.primary, subset = tissue %in% c("large intestine","rectum","small intestine"))

# print(table(adult_d10x.primary@meta.data$author_cell_type, adult_d10x.primary@meta.data$category))
# print(table(adult_d10x.primary@meta.data$author_cell_type, adult_d10x.primary@meta.data$Age_group))
# 
# print(table(adult_d10x.primary@meta.data$Age_group))

adult_d10x.primary<-subset(adult_d10x.primary, subset = Age_group %in% c("Adult"))
adult_d10x.primary


# head(colnames(adult_d10x.primary))
# head(rownames(adult_d10x.primary))

MHCI = c('ENSG00000204642', 'ENSG00000204632', 'ENSG00000206503', 'ENSG00000204592', 
         'ENSG00000204525', 'ENSG00000234745',"ENSG00000168394","ENSG00000204267",
         "ENSG00000240065","ENSG00000204264","ENSG00000166710","ENSG00000125347","ENSG00000140853")

crypt_villis = c( "ENSG00000280501" , "ENSG00000111701", "ENSG00000136826", "ENSG00000204936", "ENSG00000154269",
                  "ENSG00000084674" ,"ENSG00000196839" ,"ENSG00000137860", "ENSG00000164120" ,"ENSG00000079385", "ENSG00000171431", "ENSG00000135318",
                  "ENSG00000165949" ,"ENSG00000259823" ,"ENSG00000145287" ,"ENSG00000015520" ,"ENSG00000135549" ,"ENSG00000110244", "ENSG00000118137",
                  "ENSG00000007306" ,"ENSG00000105388" ,"ENSG00000086548" ,"ENSG00000146648", "ENSG00000073737" ,"ENSG00000155366", "ENSG00000117472")


adult_d10x.primary <- AddModuleScore(
  object = adult_d10x.primary,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score'
)

adult_d10x.primary <- AddModuleScore(
  object = adult_d10x.primary,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)


adult_d10x.primary <-FindVariableFeatures(adult_d10x.primary, selection.method = "vst", nfeatures = 2000)
adult_d10x.primary <-ScaleData(adult_d10x.primary)
adult_d10x.primary <- RunPCA(adult_d10x.primary, npcs = 30, verbose = FALSE)
adult_d10x.primary <- RunUMAP(adult_d10x.primary, reduction = "pca", dims = 1:30)

# 
# 
# DimPlot(adult_d10x.primary, reduction = "umap", pt.size=0.25, label=T)
# ggsave(file=here("figs","adult_primary_UMAP_cluster.padult_d10x.primary_shared"), w=6,h=5)
# ggsave(file=here("figs/jpeg","adult_primary_UMAP_cluster.jpeg"), w=6,h=5)
# 
# DimPlot(adult_d10x.primary, reduction = "umap", group.by = "tissue", pt.size=0.25, label=T)
# ggsave(file=here("figs","adult_primary_epi_UMAP_tissue.padult_d10x.primary_shared"), w=8,h=5)
# ggsave(file=here("figs/jpeg","adult_primary_UMAP_tissue.jpeg"), w=8,h=5)
# 
# DimPlot(adult_d10x.primary, reduction = "umap", group.by = "cell_type", pt.size=0.25, label=T)
# ggsave(file=here("figs","adult_primary_epi_UMAP_cell_type.padult_d10x.primary_shared"), w=10,h=5)
# ggsave(file=here("figs/jpeg","adult_primary_UMAP_cell_type.jpeg"), w=10,h=5)
# 
# FeaturePlot(adult_d10x.primary, features = "MHCI_score1",reduction = "umap", min.cutoff = "q9", pt.size=1)
# ggsave(file=here("figs","adult_primary_epi_UMAP_MHCI.padult_d10x.primary_shared"), w=6,h=5)
# ggsave(file=here("figs/jpeg","adult_primary_UMAP_MHCI.jpeg"), w=6,h=5)
# 
# 
# ggplot(adult_d10x.primary@meta.data, aes(reorder(cell_type, MHCI_score1, FUN = median),MHCI_score1))+
#   geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.1,aes(fill=tissue))+
#   theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),
#                       axis.text.x = element_text(angle = 45, hjust = 1))+facet_wrap(~tissue, scales="free_x")
# ggsave(file=here("figs/adult_primary_epi_MHCI_boxplot_celltype.padult_d10x.primary_shared"), w=30, h=10)
# ggsave(file=here("figs/jpeg/adult_primary_epi_MHCI_boxplot_celltype.jpeg"), w=30, h=10)



#### data frame for plotting
umap_mat<-as.data.frame(Embeddings(object = adult_d10x.primary, reduction = "umap"))#
umap_mat$cell<-rownames(umap_mat)

meta<-adult_d10x.primary@meta.data
meta$cell<-rownames(meta)

plt_UMAP<-merge(meta, umap_mat, by="cell")


###### summary plot
table(adult_d10x.primary$donor_id, adult_d10x.primary$tissue)
adult_d10x.primary

plt<-adult_d10x.primary@meta.data[,c("MHCI_score1","crypt_villis_score1","author_cell_type","tissue")]
save(plt,plt_UMAP, file=here("data/adult_intestine","plt_data_allcompartment.RData"))

#load("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/more_single_cell/plt_data_allcompartment.RData")

epithelial<-c("BEST4 enterocyte","Stem","early enterocyte","enterocyte","enteroendocrine",
              "Goblet cell", "Paneth cell","Paneth (UC only)","Tuft","TA","Paneth",
              "BEST2+ Goblet cell","BEST4+ epithelial","Colonocyte","D cells (SST+)","EC cells (TAC1+)","EECs","Enterocyte","I cells (CCK+)","K cells (GIP+)","L cells (PYY+)",
              "M/X cells (MLN/GHRL+)","Microfold cell","N cells (NTS+)","Progenitor (NEUROG3+)","Stem cells")

immune<-c("Activated B cell","activated DC" ,"Activated T","B cell","CD4 T cell",
          "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
          "IgA plasma cell","IgG plasma cell","Macrophage" ,"mast cells","Memory B cell" ,
          "Monocyte","pDC" ,"Tfh" ,
          "Treg","Mast cell",
          "Activated CD4 T","Activated CD8 T","CD8 Tmem","CLC+ Mast cell","CX3CR1+ CD8 Tmem","ILC3","LYVE1+ Macrophage","Lymphoid DC",
          "MAIT cell","Macrophages","Memory B","NK cell","Naive B","SELL+ CD4 T","SELL+ CD8 T","TRGV2 gdT","TRGV5/7 gdT","Th1","Th17","gdT","TRGV4 gdT")

stromal_and_glial<-c("S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"Lymphatic endothelial cell","Arterial endothelial cell","pericyte","Venous endothelial cell","myofibroblast", "Glial cell",
                     "Adult Glia","Contractile pericyte (PLN+)","Fetal arterial EC","LEC1 (ACKR4+)","LEC3 (ADGRG3+)","LEC4 (STAB2+)","LEC5 (CLDN11+)","LEC6 (ADAMTS4+)",
                     "Mature arterial EC","Mature venous EC","Mesothelium (PRG4+)","Monocytes","Pericyte",
                     "Stromal 1 (ADAMDEC1+)","Stromal 1 (CCL11+)","Stromal 2 (CH25H+)", 'Stromal 2 (NPY+)',"Stromal 3 (C7+)","T reticular",
                     "Transitional Stromal 3 (C3+)","arterial capillary","cycling EC","mLN Stroma (FMO2+)","myofibroblast (RSPO2+)","venous capillary")

plt_UMAP$general_type<-NA
plt_UMAP$general_type[which(plt_UMAP$author_cell_type%in%epithelial)]<-"Epithelial"
plt_UMAP$general_type[which(plt_UMAP$author_cell_type%in%immune)]<-"Immune"
plt_UMAP$general_type[which(plt_UMAP$author_cell_type%in%stromal_and_glial)]<-"Stromal and Glial"

table(plt_UMAP[which(is.na(plt_UMAP$general_type)),]$author_cell_type)[which(table(plt_UMAP[which(is.na(plt_UMAP$general_type)),]$author_cell_type)>0)]
table(plt_UMAP$general_type)


getPalette = colorRampPalette(brewer.pal(9, "Reds"))
epi_col<-getPalette(length(unique(epithelial)))
getPalette = colorRampPalette(brewer.pal(9, "Blues"))
immune_col<-getPalette(length(unique(immune)))
getPalette = colorRampPalette(brewer.pal(9, "Greens"))
stromal_col<-getPalette(length(unique(stromal_and_glial)))

stromal_col[1]<-"#006D2C"
stromal_col[2]<-"#00441B"
stromal_col[8]<-"#41ae76"
stromal_col[9]<-"#66c2a4"
show_col(stromal_col)

epi_col[1]<-"#fc4e2a"
epi_col[2]<-"#a50f15"
show_col(epi_col)

immune_col[1]<-"#67a9cf"
immune_col[2]<-"#3690c0"
immune_col[3]<-"#4292c6"
immune_col[4]<-"#9ecae1"
immune_col[9]<-"#6baed6"
immune_col[21]<-"#9ecae1"
immune_col[22]<-"#6baed6"
show_col(immune_col)

cell_cols<-c(epi_col, immune_col, stromal_col)
names(cell_cols)<-c(epithelial, immune, stromal_and_glial)

ggplot(plt_UMAP, aes(UMAP_1, UMAP_2, color=author_cell_type)) + geom_point(size=0.25) +
  theme_classic() +th_present +
  scale_color_manual(name="Cell",values = cell_cols, drop = T, guide=F)+
  xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=5, y=10, label="Immune", color="#1C6BB0")+
  annotate("text",x=-2, y=-4, label="Epithelial", color="#BB1419")+
  annotate("text",x=-2, y=-13.5, label="Stromal and Glial", color="#41AB5D")+
  annotate("text",x=-10, y=-16, label=paste("n =", prettyNum(nrow(plt_UMAP), big.mark = ",")), color="grey40")
  
ggsave(file=here("figs","adult_primary_controls_UMAP_celltype_compartment.pdf"), w=5.5,h=5)
ggsave(file=here("figs/jpeg","adult_primary_controls_UMAP_celltype_compartment.jpeg"), w=5.5,h=5)

# 
# ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_compartment.pdf"), w=4.5,h=4)
# ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_compartment.jpeg"), w=4.5,h=4,bg="white")
# 

plt_UMAP<-plt_UMAP[(order(plt_UMAP$MHCI_score1)),]
ggplot(plt_UMAP, aes(UMAP_1, UMAP_2, color=MHCI_score1))+geom_point()+
  th_present+theme_bw()+scale_color_continuous_sequential("ag_GrnYl", rev = F,na.value = "grey80")+
  xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-10, y=-16, label=paste("n =", prettyNum(nrow(plt), big.mark = ",")), color="grey40")+
  theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))
ggsave(file=here("figs", "adult_primary_epi_MHCI_by_celltype.pdf"), w=9,h=8)
ggsave(file=here("figs/jpeg", "adult_primary_epi_MHCI_celltype.jpeg"), w=9,h=8,bg="white")
# 
# ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_compartment_MHCI_Score.pdf"), w=5,h=4)
# ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_compartment_MHCI_Score.jpeg"), w=5,h=4,bg="white")
# 



### boxplots for final figure
plt_NLRC5<-plt_exp[which(plt_exp$gene_name=="NLRC5"),]
comp_simple<-list(c("Epithelial","Stromal and Glial"),c("Epithelial","Immune"),c("Stromal and Glial","Immune"))

          # ggplot(plt_NLRC5, aes(tissue, value))+geom_violin(fill="grey70", color="grey70")+
          #   geom_boxplot(aes(fill=tissue),width=0.2,outlier.shape = NA)+th_present+theme_classic()+
          #   scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), guide="none")+ylab("NLRC5 Count")+xlab("")+
          #   geom_signif(comparisons = comp_simple, step_increase = 0.08,tip_length = 0.01,
          #               size = 0.3,vjust = 0.5,textsize = 6,  map_signif_level = T, color="grey60")
          # 
          # ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_epi_NLRC5_boxplot.pdf"), w=4,h=3)
          # ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_epi_NLRC5_boxplot.jpeg"), w=5,h=4,bg="white")

ggplot(plt_UMAP, aes(general_type, MHCI_score1))+geom_violin(fill="grey70", color="grey70")+
  geom_boxplot(aes(fill=general_type),width=0.2,outlier.shape = NA)+th_present+theme_classic()+
  scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), guide="none")+ylab("MHC-I Score")+xlab("")+
  geom_signif(comparisons = comp_simple, step_increase = 0.08,tip_length = 0.01,
              size = 0.3,vjust = 0.5,textsize = 6,  map_signif_level = T, color="grey60")+facet_wrap(~tissue)
ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_epi_MHCIscore_boxplot.pdf"), w=4.5,h=2.5)
ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_epi_MHCIscore_boxplot.jpeg"), w=4.5,h=2.5,bg="white")
# 
# ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_compartment_MHCI_Score_boxplot.pdf"), w=10,h=4)
# ggsave(file=here("../EBI/MHCI/more_single_cell", "adult_primary_compartment_MHCI_Score_boxplot.jpeg"), w=10,h=4,bg="white")
# 




# 
# ### Plot UMAPS
# plt_UMAP<-plt_UMAP[which(!(plt_UMAP$author_cell_type%in%c("M/X cells (MLN/GHRL+)","Progenitor (NEUROG3+)"))),]
# plt_UMAP$author_cell_type<-as.factor(as.character(plt_UMAP$author_cell_type))
# levels(plt_UMAP$author_cell_type)<-c("BEST2+ Goblet cell","BEST4+ epithelial","Colonocyte","Enteroendocrine","Enteroendocrine",
#                                      "Enteroendocrine","Enterocyte","Goblet cell","Enteroendocrine","Enteroendocrine", "Enteroendocrine",
#                                      "Microfold cell","Enteroendocrine","Paneth","Stem cells","TA","Tuft")
# levels(plt_UMAP$author_cell_type)
# cols_manual<-c("#a30443ff", "#87d235ff","#04b026","#5c1e6cff","#208d43ff", "#ad0044ff","#5a5ce0","#d9c025ff","#6493ebff","#ef6a13ff","#a03c9fff")
# names(cols_manual) <- levels(plt_UMAP$author_cell_type)
# fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
# colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)
# 
# ggplot(plt_UMAP, aes(UMAP_1,UMAP_2, color=author_cell_type))+
#   geom_point(size=1.5)+colcale_cols_manual+
#   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))
# ggsave(file=here("figs", "adult_primary_epi_celltype.padult_d10x.primary_shared"), w=7,h=5)
# ggsave(file=here("figs/jpeg", "adult_primary_epi_celltype.jpeg"), w=7,h=5,bg="white")
# 
# ggplot(plt_UMAP, aes(UMAP_1,UMAP_2, color=tissue))+
#   geom_point(size=1.5)+scale_color_manual(values=c("#a6d96a","darkgreen","cornflowerblue"))+
#   theme_classic()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12),plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))
# ggsave(file=here("figs", "adult_primary_epi_tissue.padult_d10x.primary_shared"), w=7,h=5)
# ggsave(file=here("figs/jpeg", "adult_primary_epi_tissue.jpeg"), w=7,h=5,bg="white")
# 
# ## simplify cell type
# plt<-plt[which(!(plt$author_cell_type%in%c("M/X cells (MLN/GHRL+)","Progenitor (NEUROG3+)"))),]
# 
# ## merge all enteroendocrine types
# plt$author_cell_type<-as.factor(as.character(plt$author_cell_type))
# levels(plt$author_cell_type)<-c("BEST2+ Goblet cell","BEST4+ epithelial","Colonocyte","Enteroendocrine","Enteroendocrine",
#                                 "Enteroendocrine","Enterocyte","Goblet cell","Enteroendocrine","Enteroendocrine", "Enteroendocrine",
#                                 "Microfold cell","Enteroendocrine","Paneth","Stem cells","TA","Tuft")
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
# 
# ggsave(file=here("figs", "adult_primary_epi_MHCI_crpyt_summary.padult_d10x.primary_shared"),grid_fig, w=10,h=10)
# ggsave(file=here("figs/jpeg", "adult_primary_epi_MHCI_crpyt_summary.jpeg"),grid_fig, w=10,h=10,bg="white")
# 
# cor(plt$MHCI_score1, plt$crypt_villis_score1, method = "spearman")
# 
# 
# library(stats)
# aov(plt$MHCI_score1 ~ plt$tissue)
# 
# 
# 
# ggplot(plt, aes(tissue,MHCI_score1))+geom_violin(color="lightgrey",fill="lightgrey")+geom_boxplot(aes(fill=tissue))+facet_wrap(~author_cell_type)+scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), name="Segment")+theme_bw()+th+ylab("MHC I Score")+xlab("")
# ggsave(file=here("figs", "adult_primary_epi_MHCI_by_celltype.padult_d10x.primary_shared"), w=9,h=8)
# ggsave(file=here("figs/jpeg", "adult_primary_epi_MHCI__celltype.jpeg"), w=9,h=8,bg="white")
# 
# 
# ###############
# ## Differential expression and UMAP epithelial only individual genes
# ###############
# print("#################### DE ")
# adult_d10x.primary
# 
# d10x.exp.GOI<-FetchData(object = adult_d10x.primary, vars = MHCI)
# head(d10x.exp.GOI)
# 
# d10x.exp.GOI<-as.data.frame(t(d10x.exp.GOI))
# d10x.exp.GOI[1:5,1:5]
# d10x.exp.GOI$gene<-rownames(d10x.exp.GOI)
# d10x.exp.GOI<-melt(d10x.exp.GOI)#
# head(d10x.exp.GOI)
# 
# umap_mat<-as.data.frame(Embeddings(object = adult_d10x.primary, reduction = "umap"))#
# umap_mat$cell<-rownames(umap_mat)
# 
# meta<-adult_d10x.primary@meta.data
# meta$cell<-rownames(meta)
# 
# plt<-merge(meta, umap_mat, by="cell")
# 
# plt_exp<-merge(d10x.exp.GOI, plt,by.x="variable", by.y="cell")
# save(plt_exp, file=here("data/adult_intestine","plt_MHCIexp_data.RData"))
# #load(here("../EBI/MHCI/more_single_cell","plt_MHCIexp_data.RData"))
# 
# plt_exp$tissue<-as.factor(as.character(plt_exp$tissue))
# levels(plt_exp$tissue)
# levels(plt_exp$tissue)<-c("Large\nIntestine","Rectum","Small\nIntestine")
# 
# MHCI_df<-data.frame(ens_gene=MHCI, gene_name=c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5"))
# plt_exp<-merge(plt_exp, MHCI_df, by.x="gene", by.y="ens_gene")
# 
# 
# plt_gene_UMAP<-function(gene){
#   plt_gene<-plt_exp[which(plt_exp$gene_name==gene),]
#   print(gene)
#   print(pairwise.t.test(plt_gene$value, plt_gene$tissue,p.adjust.method = "BH", pool.sd = FALSE))
#   print(tapply(plt_gene$value,plt_gene$tissue, mean))
#   comp_simple<-list(c("Large\nIntestine","Small\nIntestine"),c("Large\nIntestine","Rectum"),c("Small\nIntestine","Rectum"))
#   
#   boxplot<-ggplot(plt_gene, aes(tissue, value))+geom_violin(fill="grey70", color="grey70")+geom_boxplot(aes(fill=tissue),width=0.2,outlier.shape = NA)+th_present+theme_bw()+
#     scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), guide="none")+ylab("Count")+xlab("Segment")+
#     geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
#                 size = 0.3,vjust = 0.5,
#                 textsize = 2,  map_signif_level = T, color="grey60")
#   
#   plt_gene<-plt_gene[(order(plt_gene$value)),]
#   UMAP_gene<-ggplot(plt_gene, aes(UMAP_1, UMAP_2, color=value))+geom_point()+
#     th_present+theme_bw()+scale_color_continuous_sequential("ag_GrnYl", rev = F,na.value = "grey80")+
#     xlab("UMAP 1")+ylab("UMAP 2")+
#     annotate("text",x=-10, y=-11, label=paste("n =", prettyNum(nrow(plt_gene), big.mark = ",")), color="grey40")
#   
#   UMAP_segemnt<-ggplot(plt_gene, aes(UMAP_1, UMAP_2, color=tissue))+geom_point()+
#     th_present+theme_bw()+scale_color_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), guide="none")+
#     xlab("UMAP 1")+ylab("UMAP 2")+
#     annotate("text",x=-10, y=-11, label=paste("n =", prettyNum(nrow(plt_gene), big.mark = ",")), color="grey40")
#   
#   plot_grid(UMAP_segemnt,UMAP_gene,boxplot, ncol=3, rel_widths = c(0.8,1,0.5))}
# 
# plt_gene_UMAP("NLRC5")
# ggsave(file=here("/home/redgar/Documents/EBI/MHCI/more_single_cell", "adult_primary_epi_MHCI_celltype_NLRC5.jpeg"), w=9,h=3,bg="white")
# 
# plt_gene_UMAP("B2M")
# ggsave(file=here("/home/redgar/Documents/EBI/MHCI/more_single_cell", "adult_primary_epi_MHCI_celltype_B2M.jpeg"), w=9,h=3,bg="white")
# 
# 
# ##save all genes
# lapply(unique(plt_exp$gene_name), function(gene){
#   plt_gene_UMAP(gene)
#   ggsave(file=here("/home/redgar/Documents/EBI/MHCI/more_single_cell", paste("adult_primary_epi_MHCI_celltype_",gene,".jpeg", sep="")), w=9,h=3,bg="white")})
# 
# 
# plt_NLRC5<-plt_exp[which(plt_exp$gene_name=="NLRC5"),]
# 
# print(pairwise.t.test(plt_NLRC5$MHCI_score1, plt_NLRC5$tissue,p.adjust.method = "BH", pool.sd = FALSE))
# print(pairwise.t.test(plt_NLRC5$value, plt_NLRC5$tissue,p.adjust.method = "BH", pool.sd = FALSE))
# 
# tapply(plt_NLRC5$value, plt_NLRC5$tissue, mean)
# 
# comp_simple<-list(c("Large\nIntestine","Small\nIntestine"),c("Large\nIntestine","Rectum"),c("Small\nIntestine","Rectum"))
# 
# plot_grid(
#   ggplot(plt_NLRC5, aes(tissue,MHCI_score1))+
#     geom_violin(fill="grey80", color="grey80")+geom_boxplot(width=0.1,aes(fill=tissue),outlier.shape = NA)+
#     theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12))+
#     scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), guide="none")+ylab("MHC I Score")+xlab("")+
#     geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
#                 size = 0.3,vjust = 0.5,
#                 textsize = 2,  map_signif_level = T, color="grey60"),
#   ggplot(plt_NLRC5, aes(tissue,value))+
#     geom_violin(fill="grey80", color="grey80")+geom_boxplot(width=0.1,aes(fill=tissue),outlier.shape = NA)+
#     theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12))+
#     scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), guide="none")+ylab("NLRC5 Count")+xlab("")+
#     geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
#                 size = 0.3,vjust = 0.5,
#                 textsize = 2,  map_signif_level = T, color="grey60"))
# 
# ggsave(file=here("../EBI/MHCI/more_single_cell","adult_intestine_boxplot_compartment_MHCI_both.pdf"), w=5,h=4)
# ggsave(file=here("../EBI/MHCI/more_single_cell","adult_intestine_boxplot_compartment_MHCI_both.jpeg"), w=5,h=4)
# 
# 
# ## simplify cell type
# plt_NLRC5<-plt_NLRC5[which(!(plt_NLRC5$author_cell_type%in%c("M/X cells (MLN/GHRL+)","Progenitor (NEUROG3+)"))),]
# 
# ## merge all enteroendocrine types
# plt_NLRC5$author_cell_type<-as.factor(as.character(plt_NLRC5$author_cell_type))
# levels(plt_NLRC5$author_cell_type)
# levels(plt_NLRC5$author_cell_type)<-c("BEST2+ Goblet cell","BEST4+ epithelial","Colonocyte","Enteroendocrine","Enteroendocrine",
#                                       "Enteroendocrine","Enterocyte","Goblet cell","Enteroendocrine","Enteroendocrine", "Enteroendocrine",
#                                       "Microfold cell","Enteroendocrine","Paneth","Stem cells","TA","Tuft")
# 
# ggplot(plt_NLRC5, aes(tissue,MHCI_score1))+
#   geom_violin(fill="grey80", color="grey80")+geom_boxplot(width=0.1,aes(fill=tissue),outlier.shape = NA)+
#   theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12))+
#   scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), guide="none")+ylab("MHC I Score")+xlab("")+
#   geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
#               size = 0.3,vjust = 0.5,
#               textsize = 2,  map_signif_level = T, color="grey60")+
#   facet_wrap(~author_cell_type)
# ggsave(file=here("../EBI/MHCI/more_single_cell","adult_intestine_boxplot_compartment_MHCI_both_cell_type.pdf"), w=7,h=10)
# ggsave(file=here("../EBI/MHCI/more_single_cell","adult_intestine_boxplot_compartment_MHCI_both_cell_type.jpeg"), w=7,h=10)
# 
# 
# ggplot(plt_NLRC5, aes(tissue,value))+
#   geom_violin(fill="grey80", color="grey80")+geom_boxplot(width=0.1,aes(fill=tissue),outlier.shape = NA)+
#   theme_bw()+th+theme(legend.text=element_text(size=10),legend.title=element_text(size=12))+
#   scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), guide="none")+ylab("NLRC5 Count")+xlab("")+
#   geom_signif(comparisons = comp_simple, step_increase = 0.03,tip_length = 0.01,
#               size = 0.3,vjust = 0.5,
#               textsize = 2,  map_signif_level = T, color="grey60")+
#   facet_wrap(~author_cell_type)
# 
# 
# #########
# ## DE
# #########
# ## simplify cell type
# adult_d10x.primary@meta.data<-adult_d10x.primary@meta.data[which(!(adult_d10x.primary@meta.data$author_cell_type%in%c("M/X cells (MLN/GHRL+)","Progenitor (NEUROG3+)"))),]
# 
# ## merge all enteroendocrine types
# adult_d10x.primary@meta.data$author_cell_type<-as.factor(as.character(adult_d10x.primary@meta.data$author_cell_type))
# levels(adult_d10x.primary@meta.data$author_cell_type)
# levels(adult_d10x.primary@meta.data$author_cell_type)<-c("BEST2+ Goblet cell","BEST4+ epithelial","Colonocyte","Enteroendocrine","Enteroendocrine",
#                                                      "Enteroendocrine","Enterocyte","Goblet cell","Enteroendocrine","Enteroendocrine", "Enteroendocrine",
#                                                      "Microfold cell","Enteroendocrine","Paneth","Stem cells","TA","Tuft")
# adult_d10x.primary@meta.data$author_cell_type<-as.character(adult_d10x.primary@meta.data$author_cell_type)
# 
# ## only shared cell types
# adult_d10x.primary_shared<-subset(adult_d10x.primary, subset = author_cell_type %in% c("BEST2+ Goblet cell","BEST4+ epithelial","Enteroendocrine","Goblet cell","Microfold cell","Paneth","Stem cells","TA","Tuft"))
# adult_d10x.primary_shared
# 
# ## testing factor
# adult_d10x.primary_shared$cell_tissue<-paste(adult_d10x.primary_shared$author_cell_type, adult_d10x.primary_shared$tissue, sep = "_")
# Idents(adult_d10x.primary_shared) <- "cell_tissue"
# 
# 
# #MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
# #consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
# #and linear regression for the continuous process (i.e., the expression level). 
# 
# cell_types<-unique(adult_d10x.primary_shared$author_cell_type)
# 
# contrasts_celltype_tissue<-rbind(
#   combinations(n = 3, r = 2, v = unique(adult_d10x.primary_shared$cell_tissue)[grep("BEST2+",unique(adult_d10x.primary_shared$cell_tissue))], repeats.allowed = FALSE),
#   combinations(n = 3, r = 2, v = unique(adult_d10x.primary_shared$cell_tissue)[grep("BEST4+",unique(adult_d10x.primary_shared$cell_tissue))], repeats.allowed = FALSE),
#   combinations(n = 3, r = 2, v = unique(adult_d10x.primary_shared$cell_tissue)[grep("Goblet",unique(adult_d10x.primary_shared$cell_tissue))], repeats.allowed = FALSE),
#   combinations(n = 3, r = 2, v = unique(adult_d10x.primary_shared$cell_tissue)[grep("Microfold",unique(adult_d10x.primary_shared$cell_tissue))], repeats.allowed = FALSE),
#   combinations(n = 3, r = 2, v = unique(adult_d10x.primary_shared$cell_tissue)[grep("TA",unique(adult_d10x.primary_shared$cell_tissue))], repeats.allowed = FALSE),
#   combinations(n = 2, r = 2, v = c("Paneth_small intestine", "Paneth_large intestine"), repeats.allowed = FALSE),
#   combinations(n = 2, r = 2, v = c("Enteroendocrine_small intestine", "Enteroendocrine_large intestine"), repeats.allowed = FALSE),
#   combinations(n = 2, r = 2, v = c("Stem cells_small intestine", "Stem cells_large intestine"), repeats.allowed = FALSE),
#   combinations(n = 2, r = 2, v = c("Tuft_small intestine", "Tuft_large intestine"), repeats.allowed = FALSE))
# 
# contrasts_celltype_tissue
# 
# nrow(contrasts_celltype_tissue)
# 
# diff_exp_all<-lapply(1:nrow(contrasts_celltype_tissue), function(x){
#   de<-FindMarkers(adult_d10x.primary_shared, ident.1 = contrasts_celltype_tissue[x,1], ident.2 = contrasts_celltype_tissue[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
#   print(paste(contrasts_celltype_tissue[x,1],"vs", contrasts_celltype_tissue[x,2],":", nrow(de), sep=" "))
#   de$gene<-rownames(de)
#   rownames(de)<-NULL
#   de<-de[,c(6,1:5)]
#   de$cell.1<-contrasts_celltype_tissue[x,1]
#   de$cell.2<-contrasts_celltype_tissue[x,2]
#   de})
# 
# head(diff_exp_all)
# 
# 
# save(diff_exp_all, file=here("data/adult_intestine","diff_exp_all_data.RData"))
# 
# load(here("../EBI/MHCI/more_single_cell/diff_exp_all_data.RData"))
# diff_exp_all<-do.call(rbind, diff_exp_all)
# 
# 
# ## Just used FDR in primary so will here too
# #sig_diff_exp<-diff_exp_all[which(diff_exp_all$p_val_adj<0.005 & abs(diff_exp_all$avg_log2FC)>1),]
# sig_diff_exp<-diff_exp_all[which(diff_exp_all$p_val_adj<0.005),]
# 
# sig_diff_exp_MHC<-merge(sig_diff_exp, MHCI_df, by.x="gene", by.y="ens_gene")
# sig_diff_exp_MHC
# unique(sig_diff_exp_MHC$gene_name)
# 
# 
# check_FC<-function(gene){
#   plt_NLRC5<-plt_exp[which(plt_exp$gene_name==gene),]
#   print(pairwise.t.test(plt_NLRC5$MHCI_score1, plt_NLRC5$tissue,p.adjust.method = "BH", pool.sd = FALSE))
#   print(pairwise.t.test(plt_NLRC5$value, plt_NLRC5$tissue,p.adjust.method = "BH", pool.sd = FALSE))
#   tapply(plt_NLRC5$value, list(plt_NLRC5$author_cell_type,plt_NLRC5$tissue), mean)}
# 
# check_FC("NLRC5")
# check_FC("B2M")
# check_FC("IRF1")
# 
# 
# ## higher expressed in large
# sig_diff_exp_MHC_up<-sig_diff_exp_MHC[which(sig_diff_exp_MHC$avg_log2FC>0),]
# sig_diff_exp_MHC_up[grep("large",sig_diff_exp_MHC_up$cell.1),]
# 
# ## lower expressed in large
# sig_diff_exp_MHC_down<-sig_diff_exp_MHC[which(sig_diff_exp_MHC$avg_log2F<0),]
# sig_diff_exp_MHC_down[grep("large",sig_diff_exp_MHC_down$cell.1),]
# 
# ## Up in large compared to small
# # B2M stem; 
# 
# 
# 
# #################
# ## all box plots
# #################
# ## simplify cell type
# plt_exp<-plt_exp[which(!(plt_exp$author_cell_type%in%c("M/X cells (MLN/GHRL+)","Progenitor (NEUROG3+)"))),]
# 
# ## merge all enteroendocrine types
# plt_exp$author_cell_type<-as.factor(as.character(plt_exp$author_cell_type))
# levels(plt_exp$author_cell_type)
# levels(plt_exp$author_cell_type)<-c("BEST2+ Goblet cell","BEST4+ epithelial","Colonocyte","Enteroendocrine","Enteroendocrine",
#                                     "Enteroendocrine","Enterocyte","Goblet cell","Enteroendocrine","Enteroendocrine", "Enteroendocrine",
#                                     "Microfold cell","Enteroendocrine","Paneth","Stem cells","TA","Tuft")
# plt_exp$author_cell_type<-as.character(plt_exp$author_cell_type)
# 
# ## only shared cell types
# plt_exp_shared<-plt_exp[which(plt_exp$author_cell_type %in% c("BEST2+ Goblet cell","BEST4+ epithelial","Enteroendocrine","Goblet cell","Microfold cell","Paneth","Stem cells","TA","Tuft")),]
# 
# sig_diff_exp_MHC$celltype<-sapply(1:nrow(sig_diff_exp_MHC), function(x) strsplit(sig_diff_exp_MHC$cell.1[x],"_")[[1]][1])
# sig_diff_exp_MHC$tissue<-sapply(1:nrow(sig_diff_exp_MHC), function(x) 
#   if(sig_diff_exp_MHC$avg_log2FC[x]>0){strsplit(sig_diff_exp_MHC$cell.1[x],"_")[[1]][2]}else{strsplit(sig_diff_exp_MHC$cell.2[x],"_")[[1]][2]})
# 
# sig_diff_exp_MHC_lrg_sml<-sig_diff_exp_MHC[which(sig_diff_exp_MHC$tissue!="rectum"),]
# 
# highlight<-sig_diff_exp_MHC_lrg_sml[,c("gene_name","celltype","tissue")]
# colnames(highlight)<-c("gene_name","author_cell_type","tissue")
# 
# 
# ggplot()+geom_violin(aes(tissue, value, fill=tissue),plt_exp_shared, fill="grey", color="grey")+
#   geom_boxplot(aes(tissue, value, fill=tissue),plt_exp_shared, width=0.25, outlier.shape=NA)+
#   facet_grid(gene_name~author_cell_type, scales="free_y")+theme_bw()+th_present+
#   scale_fill_manual(values=c("#a6d96a","darkgreen","cornflowerblue"), guide="none")+
#   scale_color_manual(values=c("#53d420","cornflowerblue"), guide="none")+
#   theme(legend.position = "none")+ylab("Expression Level")+xlab("Segment")+
#   geom_rect(data = highlight, 
#             fill = NA, aes(colour = tissue), xmin = -Inf,xmax = Inf,
#             ymin = -Inf,ymax = Inf, size=1.75)
# 
# 
# ggsave(file=here("/home/redgar/Documents/EBI/MHCI/more_single_cell","adult_segment_MHC_singlecell.pdf"), w=8, h=10)
# ggsave(file=here("/home/redgar/Documents/EBI/MHCI/more_single_cell","adult_segment_MHC_singlecell.jpeg"), w=8, h=10)
