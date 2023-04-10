library(RColorBrewer)

load(here("../../codon/scRNAseq_codon/data/Primary_UMAP_controls_plus_MHCI_raw.RData"))


epithelial<-c("BEST4 enterocyte","Stem","early enterocyte","enterocyte","enteroendocrine",
  "Goblet cell", "Paneth cell","Paneth (UC only)","Tuft","TA","Paneth")
immune<-c("Activated B cell","activated DC" ,"Activated T","B cell","CD4 T cell",
  "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
  "IgA plasma cell","IgG plasma cell","Macrophage" ,"mast cells","Memory B cell" ,
  "Monocyte","pDC" ,"Tfh" ,
  "Treg","Mast cell")
stromal_and_glial<-c("S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"Lymphatic endothelial cell","Arterial endothelial cell","pericyte","Venous endothelial cell","myofibroblast", "Glial cell")

plt$general_type<-NA
plt$general_type[which(plt$cluster_ID%in%epithelial)]<-"Epithelial"
plt$general_type[which(plt$cluster_ID%in%immune)]<-"Immune"
plt$general_type[which(plt$cluster_ID%in%stromal_and_glial)]<-"Stromal and Glial"


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

ggplot(plt, aes(UMAP_1, UMAP_2, color=cluster_ID)) + geom_point(size=0.25) + 
  theme_classic() +th_present +
  scale_color_manual(name="Cell",values = cell_cols, drop = T, guide=F)+
  xlab("UMAP 1")+ylab("UMAP 2")+
  annotate("text",x=-2, y=-4, label="Immune", color="#1C6BB0")+  
  annotate("text",x=10, y=-10, label="Epithelial", color="#BB1419")+
  annotate("text",x=5, y=14, label="Stromal and Glial", color="#41AB5D")

ggsave(file=here("../../codon/scRNAseq_codon/figs","primary_controls_UMAP_celltype_compartment.pdf"), w=5.5,h=5)
ggsave(file=here("../../codon/scRNAseq_codon/figs/jpeg","primary_controls_UMAP_celltype_compartment.jpeg"), w=5.5,h=5)


#' 
#' 
#' levels(plt$general_type)<-c(
#'   "BEST4 enterocyte","Stem","early enterocyte","enterocyte","enteroendocrine",
#'   "Goblet cell", "Paneth cell","Paneth (UC only)","Tuft",
#'   "Activated B cell","activated DC" ,"Activated T","Arterial endothelial cell" ,"B cell","CD4 T cell",
#'   "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
#'   "IgA plasma cell","IgG plasma cell","Lymphatic endothelial cell","Macrophage" ,"mast cells","Memory B cell" ,
#'   "Monocyte","pDC" ,"pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"TA","Tfh" ,
#'   "Treg", "Venous endothelial cell","Glial cell","myofibroblast","Mast cell")
#' 
#' #' ### dev cell labels
#' cols_manual<-c(
#'   "#3D0AF2","#6521C1","#67D364","#367C34",
#'   "#B921C1","#308AC8",
#'   "#C86E30","#CFA218","#C12134","#238443","#810f7c" ,"#02818a","#8c510a" ,"#78c679","#67a9cf",
#'   "#3690c0","#8073ac","#b2abd2","#238443","#dd3497","#1d91c0","#addd8e","#014636",
#'   "#253494","#081d58","#bf812d","#7a0177" ,"#ce1256","#006d2c" ,
#'   "#a50f15","#542788" ,"#e31a1c","#e08214","#ef6548","#fdd49e" ,"#f46d43","#4575b4" ,
#'   "#016c59", "#bf812d","#c51b7d","#fb9a99","#fb6a4a")
