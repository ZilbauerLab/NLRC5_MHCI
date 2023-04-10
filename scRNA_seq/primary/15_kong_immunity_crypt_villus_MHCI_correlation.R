### Load libraries
library(here)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scales)
library(gridExtra)
library(reshape2)
library(gtools)
library(colorspace)
library(cowplot)
library(ggsignif)



#source(here("../EBI/scRNAseq_codon/scripts/00_pretty_plots.R"))


# dataset_loc <- here("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/expression")
dataset_loc <- here("data/kong_immunity/SCP1884/expression")

samples<-list.files(dataset_loc,recursive=T)
print(samples)
sample_IDs<-gsub(".scp.barcodes.tsv","",sapply(1:length(samples[grep("barcodes.tsv",samples)]), function(x)  strsplit(samples[grep("barcodes.tsv",samples)][x],"/")[[1]][2]))
sample_IDs

samples_notcopy<-samples[grep("barcodes.tsv",samples)]
samples_notcopy[grep("CO_STR|CO_EPI|CO_IMM|TI_STR|TI_EPI|TI_IMM",samples_notcopy)]
folders<-sapply(1:length(samples_notcopy), function(x)  strsplit(samples_notcopy[x],"/")[[1]][1])

#' 
#' 
#' d10x.list <- sapply(1:length(sample_IDs), function(y){
#'   print(sample_IDs[y])
#'   print(file.path(dataset_loc,folders[y]))
#'   mtx <- paste(file.path(dataset_loc,folders[y]),"/",sample_IDs[y],".scp.raw.mtx", sep="")
#'   cells <- paste(file.path(dataset_loc,folders[y]),"/",sample_IDs[y],".scp.barcodes.tsv", sep="")
#'   features <- paste(file.path(dataset_loc,folders[y]),"/",sample_IDs[y],".scp.features.tsv", sep="")
#'   counts <- ReadMtx(mtx = mtx, cells = cells, features = features)
#'   
#'   #' Initialize the Seurat object with the raw (non-normalized data).
#'   d10x<-CreateSeuratObject(counts = counts, project = sample_IDs[y], min.cells = 0, min.features = 0)
#'   d10x$orig.ident<-sample_IDs[y]
#'   d10x
#' })
#' 
#' d10x.list
#' 
#' d10x <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "adult_ped_map")#add.cell.ids = alldata_names2,
#' 
#' d10x
#' 
#' saveRDS(d10x, file = here("data/kong_immunity","d10x_kong_immunity.rds"))
#' 
#' 
#' ###############
#' ## only epithelial
#' ##############
#' sample_IDs_epi<-sample_IDs[grep("EPI",sample_IDs)]
#' samples_notcopy_epi<-samples_notcopy[grep("EPI",samples_notcopy)]
#' folders_epi<-sapply(1:length(samples_notcopy_epi), function(x)  strsplit(samples_notcopy_epi[x],"/")[[1]][1])
#' 
#' d10x.list <- sapply(1:length(sample_IDs_epi), function(y){
#'   print(sample_IDs_epi[y])
#'   print(file.path(dataset_loc,folders_epi[y]))
#'   mtx <- paste(file.path(dataset_loc,folders_epi[y]),"/",sample_IDs_epi[y],".scp.raw.mtx", sep="")
#'   cells <- paste(file.path(dataset_loc,folders_epi[y]),"/",sample_IDs_epi[y],".scp.barcodes.tsv", sep="")
#'   features <- paste(file.path(dataset_loc,folders_epi[y]),"/",sample_IDs_epi[y],".scp.features.tsv", sep="")
#'   counts <- ReadMtx(mtx = mtx, cells = cells, features = features)
#'   
#'   #' Initialize the Seurat object with the raw (non-normalized data).
#'   d10x<-CreateSeuratObject(counts = counts, project = sample_IDs_epi[y], min.cells = 0, min.features = 0)
#'   d10x$orig.ident<-sample_IDs_epi[y]
#'   d10x
#' })
#' 
#' d10x.list
#' 
#' d10x_epi <- merge(d10x.list[[1]], y= d10x.list[2:length(d10x.list)], merge.data=TRUE, project = "adult_ped_map")#add.cell.ids = alldata_names2,
#' 
#' d10x_epi
#' 
#' saveRDS(d10x_epi, file = here("data/kong_immunity","d10x_kong_immunity_epithelial.rds"))
#' 
#' 

####################
## Subset to just the samples they annotated 
####################
## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x<-readRDS(here("data/kong_immunity","d10x_kong_immunity.rds"))

#meta<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/metadata/scp_metadata_combined.v2.txt",header=T)
meta<-read.table(here("data/kong_immunity/SCP1884/metadata/scp_metadata_combined.v2.txt"),header=T)
meta <- meta[-(1),]

meta<-meta[which(meta$NAME%in%rownames(d10x@meta.data)),]
meta<-meta[match(meta$NAME, rownames(d10x@meta.data)),]
identical(meta$NAME, rownames(d10x@meta.data))

rownames(meta)<-meta$NAME

## add meta
d10x <- AddMetaData(d10x, metadata = meta)


####################
## Score correlation plot 
####################
MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5")
crypt_villis = c("SEPP1", "CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8", "APOBEC1", "APOB", "APOA4", "APOA1", "NPC1L1", "EGFR", "KLF4", "ENPP3", "NT5E", "SLC28A2", "ADA")


##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x <- NormalizeData(d10x,scale.factor = 10000, normalization.method = "LogNormalize")
d10x_exp<-FetchData(object = d10x, vars = c(MHCI,crypt_villis))


d10x <- AddModuleScore(
  object = d10x,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score'
)

d10x <- AddModuleScore(
  object = d10x,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)

score_data<-d10x@meta.data[,c("MHCI_score1","crypt_villis_score1","NAME")]
score_data_allcompartment<-score_data
colnames(score_data_allcompartment)<-c("MHCI_score1_allcompartment","crypt_villis_score1_allcompartment","NAME")


################
## Just in Epithelial
################
d10x_epi<-readRDS(here("data/kong_immunity","d10x_kong_immunity_epithelial.rds"))

#meta<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/metadata/scp_metadata_combined.v2.txt",header=T)
meta<-read.table(here("data/kong_immunity/SCP1884/metadata/scp_metadata_combined.v2.txt"),header=T)
meta <- meta[-(1),]

meta<-meta[which(meta$NAME%in%rownames(d10x_epi@meta.data)),]
meta<-meta[match(meta$NAME, rownames(d10x_epi@meta.data)),]
identical(meta$NAME, rownames(d10x_epi@meta.data))

rownames(meta)<-meta$NAME


## add meta
d10x_epi <- AddMetaData(d10x_epi, metadata = meta)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x_epi <- NormalizeData(d10x_epi,scale.factor = 10000, normalization.method = "LogNormalize")
d10x_exp_epi<-FetchData(object = d10x_epi, vars = c(MHCI,crypt_villis))


d10x_epi <- AddModuleScore(
  object = d10x_epi,
  features = list(MHCI),
  ctrl = 5,
  name = 'MHCI_score'
)

d10x_epi <- AddModuleScore(
  object = d10x_epi,
  features = list(crypt_villis),
  ctrl = 5,
  name = 'crypt_villis_score'
)

score_data<-d10x_epi@meta.data[,c("MHCI_score1","crypt_villis_score1","NAME")]
score_data_epithelial<-score_data
colnames(score_data_epithelial)<-c("MHCI_score1_epithelial","crypt_villis_score1_epithelial","NAME")


save(d10x_exp_epi, score_data_epithelial, d10x_exp, score_data_allcompartment, file=here("data/kong_immunity","Kong_Immunity_epithelial_MHCI_raw.RData"))

#load(here("../EBI/scRNAseq_codon/data/kong_immunity","Kong_Immunity_epithelial_MHCI_raw.RData"))


###############
## DGE
###############
print("DGE analysis")
d10x_epi$disease__ontology_label<-as.factor(d10x_epi$disease__ontology_label)
levels(d10x_epi$disease__ontology_label)<-c("CD","Normal")

#######
## split sites
#######
d10x_epi_CO<-subset(d10x_epi, subset = Site == "CO")
d10x_epi_TI<-subset(d10x_epi, subset = Site == "TI")

d10x_epi_TI$Celltype<-as.factor(d10x_epi_TI$Celltype)

d10x_epi_TI<-subset(d10x_epi_TI, subset = Celltype != "Epithelial HBB HBA")
d10x_epi_TI$Celltype<-as.factor(as.character(d10x_epi_TI$Celltype))
print(levels(d10x_epi_TI$Celltype))
levels(d10x_epi_TI$Celltype)<-c("Enteroendocrine","BEST4 Enterocyte","Enterocyte","Enterocyte","Enterocyte", 
                             "TA","Goblet" ,"Goblet" , "Goblet","Enteroendocrine",
                             "Paneth","Stem"  , "Stem",   
                             "Stem","Stem","Tuft")

d10x_epi_CO$Celltype<-as.factor(d10x_epi_CO$Celltype)
d10x_epi_CO$Celltype<-as.factor(as.character(d10x_epi_CO$Celltype))
print(levels(d10x_epi_CO$Celltype))
levels(d10x_epi_CO$Celltype)<-c("BEST4 Enterocyte","Enterocyte",
                                "Enterocyte", "Enteroendocrine",
                                "TA","Goblet" ,"Goblet" , "Goblet",
                                "Paneth","Stem"  ,      
                                "Stem","Stem","Tuft")



## testing factor TI
d10x_epi_TI$cell_stim<-paste(d10x_epi_TI$Celltype, d10x_epi_TI$disease__ontology_label, sep = "_")
Idents(d10x_epi_TI) <- "cell_stim"
table(d10x_epi_TI$Celltype, d10x_epi_TI$disease__ontology_label)

#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 
cell_types<-unique(d10x_epi_TI$Celltype)
contrasts_celltype_TI<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x_epi_TI$cell_stim[grep(cell_types[x],d10x_epi_TI$cell_stim)], repeats.allowed = FALSE)}))
contrasts_celltype_TI
nrow(contrasts_celltype_TI)


## testing factor CO
d10x_epi_CO$cell_stim<-paste(d10x_epi_CO$Celltype, d10x_epi_CO$disease__ontology_label, sep = "_")
Idents(d10x_epi_CO) <- "cell_stim"
table(d10x_epi_CO$Celltype, d10x_epi_CO$disease__ontology_label)

#MAST (Finak et al., 2015), which fits a hurdle model to the expression of each gene,
#consisting of logistic regression for the zero process (i.e., whether the gene is expressed) #
#and linear regression for the continuous process (i.e., the expression level). 
cell_types<-unique(d10x_epi_CO$Celltype)
contrasts_celltype_CO<-do.call(rbind,lapply(1:length(cell_types), function(x){
  combinations(n = 2, r = 2, v = d10x_epi_CO$cell_stim[grep(cell_types[x],d10x_epi_CO$cell_stim)], repeats.allowed = FALSE)}))
contrasts_celltype_CO
nrow(contrasts_celltype_CO)


diff_exp_CO<-lapply(1:nrow(contrasts_celltype_TI), function(x){
  de<-FindMarkers(d10x_epi_CO, ident.1 = contrasts_celltype_TI[x,1], ident.2 = contrasts_celltype_TI[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  print(paste(contrasts_celltype_TI[x,1],"vs", contrasts_celltype_TI[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_TI[x,1]
  de$cell.2<-contrasts_celltype_TI[x,2]
  de})

diff_exp_CO<-do.call(rbind, diff_exp_CO)

diff_exp_TI<-lapply(1:nrow(contrasts_celltype_CO), function(x){
  de<-FindMarkers(d10x_epi_TI, ident.1 = contrasts_celltype_CO[x,1], ident.2 = contrasts_celltype_CO[x,2], test.use = "MAST",latent.vars="nFeature_RNA", verbose=F)
  print(paste(contrasts_celltype_CO[x,1],"vs", contrasts_celltype_CO[x,2],":", nrow(de), sep=" "))
  de$gene<-rownames(de)
  rownames(de)<-NULL
  de<-de[,c(6,1:5)]
  de$cell.1<-contrasts_celltype_CO[x,1]
  de$cell.2<-contrasts_celltype_CO[x,2]
  de})

diff_exp_TI<-do.call(rbind, diff_exp_TI)

save(diff_exp_CO, diff_exp_TI, file=here("data/kong_immunity","kong_immunity_diff_genes.RData"))

#load(here("../EBI/scRNAseq_codon/data/kong_immunity","kong_immunity_diff_genes.RData"))

# add a significant threshold here for adjusted as not all in these lists are
diff_exp_sig_CO<-diff_exp_CO[which(diff_exp_CO$p_val_adj<0.005),]
diff_exp_sig_CO[which(diff_exp_sig_CO$gene%in%MHCI),]

diff_exp_sig_TI<-diff_exp_TI[which(diff_exp_TI$p_val_adj<0.005),]
diff_exp_sig_TI[which(diff_exp_sig_TI$gene%in%MHCI),]








######################
## Correlation plot
######################
load(here("../EBI/scRNAseq_codon/data/kong_immunity","Kong_Immunity_epithelial_MHCI_raw.RData"))

TI_UMAP<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/cluster/TI_EPI.scp.X_umap.coords.txt", header=T)
TI_UMAP<-TI_UMAP[-1,]
CO_UMAP<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/cluster/CO_EPI.scp.X_umap.coords.txt", header=T)
CO_UMAP<-CO_UMAP[-1,]

meta<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/metadata/scp_metadata_combined.v2.txt",header=T)


##############
### TI
##############
plt<-merge(meta, TI_UMAP, by="NAME")
plt<-merge(plt, score_data_allcompartment, by="NAME")

plt<-plt[which(plt$Site == "TI"),]

            # how many samples
            length(unique(plt$donor_id))
            plt_count<-plt[,c("donor_id","disease__ontology_label")]
            plt_count<-plt_count[!duplicated(plt_count),]
            table(plt_count$disease__ontology_label)
            tapply(plt_count$donor_id, plt_count$disease__ontology_label, function(x) length(unique(x)))

plt<-plt[which(plt$Celltype != "Epithelial HBB HBA"),]
plt$Celltype<-as.factor(plt$Celltype)
levels(plt$Celltype)<-c("Enteroendocrine","BEST4 Enterocyte","Enterocyte","Enterocyte","Enterocyte",
                                "TA","Goblet" ,"Goblet" , "Goblet","Enteroendocrine",
                                "Paneth","Stem"  , "Stem",
                                "Stem","Stem","Tuft")




cols_manual_less<-c( "#5e1e6fff",
                     "#87d435ff",
                     "#238b43ff","#f16916ff","#aa0044ff",
                     "#d9c026ff","cornflowerblue",
                     "#a13da1ff")

plt$disease__ontology_label<-as.factor(plt$disease__ontology_label)
levels(plt$disease__ontology_label)<-c("CD","Control")

plt<-plt[order(rev(plt$X)),]
scatter<-ggplot(plt, aes(crypt_villis_score1_allcompartment, MHCI_score1_allcompartment, fill=disease__ontology_label, color=disease__ontology_label))+geom_point(shape=21, color="black", size=1.5)+
  scale_fill_manual(values=c("dodgerblue3","lightgrey"))+scale_color_manual(values=c("dodgerblue","grey80"))+
  theme_bw()+th_present+theme(legend.position = "none",plot.margin = unit(c(0,0,3.95,4.2), "cm"))+
  stat_smooth(method="lm",se=F)+xlab("Crypt-Villus Score")+ylab("MHC I Score")+
  annotate("text", x=1.1, y=-0.75, label=paste("Rs = ",signif(cor.test(plt$crypt_villis_score1, plt$MHCI_score1, method="spearman")$estimate,2),";",
                                               "p = ",signif(cor.test(plt$crypt_villis_score1, plt$MHCI_score1, method="spearman")$p.value,1)))



violin_celltype<-ggplot(plt, aes(reorder(Celltype, crypt_villis_score1_allcompartment), crypt_villis_score1_allcompartment, fill=Celltype))+
  geom_violin()+scale_x_discrete(position = "bottom") +coord_flip()+
  scale_fill_manual(values=cols_manual_less)+
  theme_bw()+th_present+ylab("")+xlab("")+theme(axis.title.x=element_blank(),
                                                axis.text.x=element_blank(),
                                                axis.ticks.x=element_blank(),
                                                legend.position = "none",
                                                plot.margin = unit(c(1.1,0,0,1), "cm"))

comp_simple<-list(c("CD","Control"))

box_diagnosis<-ggplot(plt, aes(disease__ontology_label, MHCI_score1_allcompartment, fill=disease__ontology_label))+geom_boxplot()+
  scale_fill_manual(values=c("dodgerblue3","lightgrey"), name="Diagnosis")+xlab("Diagnosis")+
  theme_bw()+th_present+theme(axis.title.y =element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank(),
                              plot.margin = unit(c(0,0,3,0), "cm"),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_signif(comparisons = comp_simple, step_increase = 0.05,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")

mn_cryptvillus<-tapply(plt$crypt_villis_score1_allcompartment, plt$Celltype, mean)
cryptvillus_order<-names(mn_cryptvillus)[(order(mn_cryptvillus))]

plt_count<-melt(table(plt$disease__ontology_label, plt$Celltype))
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

ggsave(file=here("../EBI/scRNAseq_codon/figs", "MHCI_crpyt_summary_Kong_immunity_TI.pdf"),grid_fig, w=10,h=10)
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg", "MHCI_crpyt_summary_Kong_immunity_TI.jpeg"),grid_fig, w=10,h=10)

t.test(plt$MHCI_score1_allcompartment~plt$disease__ontology_label)


##############
### CO
##############
plt<-merge(meta, CO_UMAP, by="NAME")
plt<-merge(plt, score_data_allcompartment, by="NAME")

plt$Celltype<-as.factor(plt$Celltype)
levels(plt$Celltype)<-c("BEST4 Enterocyte","Enterocyte",
                        "Enterocyte", "Enteroendocrine",
                        "TA","Goblet" ,"Goblet" , "Goblet",
                        "Paneth","Stem"  ,
                        "Stem","Stem","Tuft")



cols_manual_less<-c("#87d435ff",
                     "#238b43ff","#5e1e6fff","#f16916ff","#aa0044ff",
                     "#d9c026ff","cornflowerblue",
                     "#a13da1ff")

plt$disease__ontology_label<-as.factor(plt$disease__ontology_label)
levels(plt$disease__ontology_label)<-c("CD","Control")

plt<-plt[order(rev(plt$X)),]
scatter<-ggplot(plt, aes(crypt_villis_score1_allcompartment, MHCI_score1_allcompartment, fill=disease__ontology_label, color=disease__ontology_label))+geom_point(shape=21, color="black", size=1.5)+
  scale_fill_manual(values=c("dodgerblue3","lightgrey"))+scale_color_manual(values=c("dodgerblue","grey80"))+
  theme_bw()+th_present+theme(legend.position = "none",plot.margin = unit(c(0,0,3.95,4), "cm"))+
  stat_smooth(method="lm",se=F)+xlab("Crypt-Villus Score")+ylab("MHC I Score")+
  annotate("text", x=1.1, y=-0.75, label=paste("Rs = ",signif(cor.test(plt$crypt_villis_score1, plt$MHCI_score1, method="spearman")$estimate,2),";",
                                               "p = ",signif(cor.test(plt$crypt_villis_score1, plt$MHCI_score1, method="spearman")$p.value,1)))


violin_celltype<-ggplot(plt, aes(reorder(Celltype, crypt_villis_score1_allcompartment), crypt_villis_score1_allcompartment, fill=Celltype))+
  geom_violin()+scale_x_discrete(position = "bottom") +coord_flip()+
  scale_fill_manual(values=cols_manual_less)+
  theme_bw()+th_present+ylab("")+xlab("")+theme(axis.title.x=element_blank(),
                                                axis.text.x=element_blank(),
                                                axis.ticks.x=element_blank(),
                                                legend.position = "none",
                                                plot.margin = unit(c(1.1,0,0,1), "cm"))


comp_simple<-list(c("CD","Control"))

box_diagnosis<-ggplot(plt, aes(disease__ontology_label, MHCI_score1_allcompartment, fill=disease__ontology_label))+geom_boxplot()+
  scale_fill_manual(values=c("dodgerblue3","lightgrey"), name="Diagnosis")+xlab("Diagnosis")+
  theme_bw()+th_present+theme(axis.title.y =element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank(),
                              plot.margin = unit(c(0,0,3,0), "cm"),
                              axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_signif(comparisons = comp_simple, step_increase = 0.05,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")

mn_cryptvillus<-tapply(plt$crypt_villis_score1_allcompartment, plt$Celltype, mean)
cryptvillus_order<-names(mn_cryptvillus)[(order(mn_cryptvillus))]

plt_count<-melt(table(plt$disease__ontology_label, plt$Celltype))
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

ggsave(file=here("../EBI/scRNAseq_codon/figs", "MHCI_crpyt_summary_Kong_immunity_CO.pdf"),grid_fig, w=10,h=10)
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg", "MHCI_crpyt_summary_Kong_immunity_CO.jpeg"),grid_fig, w=10,h=10)


t.test(plt$MHCI_score1_allcompartment~plt$disease__ontology_label)





#########
## MHCI boxplot and UMAP TI
#########
plt<-merge(meta, TI_UMAP, by="NAME")
plt<-merge(plt, score_data_allcompartment, by="NAME")

plt<-plt[which(plt$Site == "TI"),]

plt<-plt[which(plt$Celltype != "Epithelial HBB HBA"),]
plt$Celltype<-as.factor(plt$Celltype)
levels(plt$Celltype)<-c("Enteroendocrine","BEST4 Enterocyte","Enterocyte","Enterocyte","Enterocyte",
                        "TA","Goblet" ,"Goblet" , "Goblet","Enteroendocrine",
                        "Paneth","Stem"  , "Stem",
                        "Stem","Stem","Tuft")

## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

plt$X<-as.numeric(plt$X)
plt$Y<-as.numeric(plt$Y)

cluster_ID_plt_UMAP<-ggplot(plt, aes(X,Y, color=Celltype))+
  geom_point(size=1.5)+xlab("UMAP 1")+ylab("UMAP 2")+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  scale_color_manual(values=cols_manual_less, name="Cell Type")
leg_umap_clust = get_leg(cluster_ID_plt_UMAP)

cluster_ID_plt_UMAP<- grid.arrange(cluster_ID_plt_UMAP+theme(legend.position = "none"),
                                   leg_umap_clust, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_immunity_aligned_cell_label_UMAP_TI.jpeg"),cluster_ID_plt_UMAP, w=9, h=7.5)

plt<-plt[order(plt$n_counts),]
stim_plt_UMAP<-ggplot(plt, aes(X,Y, color=disease__ontology_label))+geom_point(size=1)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  scale_color_manual(values=c("dodgerblue3","lightgrey"), name="Diagnosis")+xlab("UMAP 1")+ylab("UMAP 2")
leg_umap_stim = get_leg(stim_plt_UMAP)

stim_plt_UMAP<- grid.arrange(stim_plt_UMAP+theme(legend.position = "none"),
                             leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_immunity_aligned_diagnosis_UMAP_TI.jpeg"),stim_plt_UMAP,w=9, h=7.5)


MHC_plt_UMAP<-ggplot(plt, aes(X,Y, color=MHCI_score1_allcompartment))+geom_point(size=1)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  labs(color = "Pathway\nActivation\nScore")+xlab("UMAP 1")+ylab("UMAP 2")
leg_umap_stim = get_leg(MHC_plt_UMAP)

MHC_plt_UMAP<- grid.arrange(MHC_plt_UMAP+theme(legend.position = "none"),
                            leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_immunity_aligned_MHC1_UMAP_TI.jpeg"),MHC_plt_UMAP, w=9, h=7.5)



crypt_villis_plt_UMAP<-ggplot(plt, aes(X,Y, color=crypt_villis_score1_allcompartment))+geom_point(size=1)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  labs(color = "Pathway\nActivation\nScore")+xlab("UMAP 1")+ylab("UMAP 2")
leg_umap_stim = get_leg(crypt_villis_plt_UMAP)

crypt_villis_plt_UMAP<- grid.arrange(crypt_villis_plt_UMAP+theme(legend.position = "none"),
                                     leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_immunity_aligned_crypt_villis_UMAP_TI.jpeg"),crypt_villis_plt_UMAP,w=9, h=7.5)



comp_simple<-list(c("Crohn's disease","normal"))
plt$Celltype<-factor(plt$Celltype, levels=rev(c("Stem","Paneth","TA","Tuft","Enteroendocrine","BEST4 Enterocyte","Enterocyte", "Goblet" )))

ggplot(plt, aes(disease__ontology_label,MHCI_score1_allcompartment))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.2,aes(fill=disease__ontology_label))+
  theme_bw()+th+
  facet_wrap(~Celltype, ncol=2)+  scale_fill_manual(values=c("cornflowerblue","lightgrey"), name="Treatment")+
  ylab("MHC I Score")+xlab("")+
  geom_signif(comparisons = comp_simple, step_increase = 0.05,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")+
  th_present+theme(strip.background = element_rect(fill = "white"),
                   legend.position = "none")
ggsave(file=here("../EBI/scRNAseq_codon/figs","Kong_immunity_MHCI_score_boxplot_TI.pdf"),  w=4, h=7)
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_immunity_MHCI_score_boxplot_TI.jpeg"),  w=4, h=7)

tapply(plt$MHCI_score1_allcompartment, list(plt$Celltype, plt$disease__ontology_label), mean)



#########
## MHCI boxplot and UMAP CO
#########
plt<-merge(meta, CO_UMAP, by="NAME")
plt<-merge(plt, score_data_allcompartment, by="NAME")

plt$Celltype<-as.factor(plt$Celltype)
levels(plt$Celltype)<-c("BEST4 Enterocyte","Enterocyte",
                        "Enterocyte", "Enteroendocrine",
                        "TA","Goblet" ,"Goblet" , "Goblet",
                        "Paneth","Stem"  ,
                        "Stem","Stem","Tuft")

## grab legened from plot
get_leg = function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}

plt$X<-as.numeric(plt$X)
plt$Y<-as.numeric(plt$Y)

cluster_ID_plt_UMAP<-ggplot(plt, aes(X,Y, color=Celltype))+
  geom_point(size=1.5)+xlab("UMAP 1")+ylab("UMAP 2")+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  scale_color_manual(values=cols_manual_less, name="Cell Type")
leg_umap_clust = get_leg(cluster_ID_plt_UMAP)

cluster_ID_plt_UMAP<- grid.arrange(cluster_ID_plt_UMAP+theme(legend.position = "none"),
                                   leg_umap_clust, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_immunity_aligned_cell_label_UMAP_CO.jpeg"),cluster_ID_plt_UMAP, w=9, h=7.5)

plt<-plt[order(plt$n_counts),]
stim_plt_UMAP<-ggplot(plt, aes(X,Y, color=disease__ontology_label))+geom_point(size=1)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  scale_color_manual(values=c("dodgerblue3","lightgrey"), name="Diagnosis")+xlab("UMAP 1")+ylab("UMAP 2")
leg_umap_stim = get_leg(stim_plt_UMAP)

stim_plt_UMAP<- grid.arrange(stim_plt_UMAP+theme(legend.position = "none"),
                             leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_immunity_aligned_diagnosis_UMAP_CO.jpeg"),stim_plt_UMAP,w=9, h=7.5)


MHC_plt_UMAP<-ggplot(plt, aes(X,Y, color=MHCI_score1_allcompartment))+geom_point(size=1)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  labs(color = "Pathway\nActivation\nScore")+xlab("UMAP 1")+ylab("UMAP 2")
leg_umap_stim = get_leg(MHC_plt_UMAP)

MHC_plt_UMAP<- grid.arrange(MHC_plt_UMAP+theme(legend.position = "none"),
                            leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_immunity_aligned_MHC1_UMAP_CO.jpeg"),MHC_plt_UMAP, w=9, h=7.5)



crypt_villis_plt_UMAP<-ggplot(plt, aes(X,Y, color=crypt_villis_score1_allcompartment))+geom_point(size=1)+
  theme_classic()+th_present+theme(plot.margin = unit(c(0.5,0,0.5,0.7), "cm"))+
  labs(color = "Pathway\nActivation\nScore")+xlab("UMAP 1")+ylab("UMAP 2")
leg_umap_stim = get_leg(crypt_villis_plt_UMAP)

crypt_villis_plt_UMAP<- grid.arrange(crypt_villis_plt_UMAP+theme(legend.position = "none"),
                                     leg_umap_stim, ncol=2, widths=c(0.8,0.2))
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_immunity_aligned_crypt_villis_UMAP_CO.jpeg"),crypt_villis_plt_UMAP,w=9, h=7.5)



comp_simple<-list(c("Crohn's disease","normal"))
plt$Celltype<-factor(plt$Celltype, levels=rev(c("Stem","Paneth","TA","Tuft","Enteroendocrine","BEST4 Enterocyte","Enterocyte", "Goblet" )))

ggplot(plt, aes(disease__ontology_label,MHCI_score1_allcompartment))+
  geom_violin(fill="grey80", color="white")+geom_boxplot(width=0.2,aes(fill=disease__ontology_label))+
  theme_bw()+th+
  facet_wrap(~Celltype, ncol=2)+  scale_fill_manual(values=c("cornflowerblue","lightgrey"), name="Treatment")+
  ylab("MHC I Score")+xlab("")+
  geom_signif(comparisons = comp_simple, step_increase = 0.05,tip_length = 0.01,
              size = 0.3,vjust = 0.5,
              textsize = 3,  map_signif_level = T, color="grey60")+
  th_present+theme(strip.background = element_rect(fill = "white"),
                   legend.position = "none")
ggsave(file=here("../EBI/scRNAseq_codon/figs","Kong_immunity_MHCI_score_boxplot_CO.pdf"),  w=4, h=7)
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_immunity_MHCI_score_boxplot_CO.jpeg"),  w=4, h=7)


tapply(plt$MHCI_score1_allcompartment, list(plt$Celltype, plt$disease__ontology_label), mean)







############
## Heat plot MHCI TI
############
load(here("../EBI/scRNAseq_codon/data/kong_immunity","Kong_Immunity_epithelial_MHCI_raw.RData"))
TI_UMAP<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/cluster/TI_EPI.scp.X_umap.coords.txt", header=T)
TI_UMAP<-TI_UMAP[-1,]
meta<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/metadata/scp_metadata_combined.v2.txt",header=T)

plt<-merge(meta, TI_UMAP, by="NAME")

plt<-plt[which(plt$Celltype != "Epithelial HBB HBA"),]
plt$Celltype<-as.factor(plt$Celltype)
levels(plt$Celltype)<-c("Enteroendocrine","BEST4 Enterocyte","Enterocyte","Enterocyte","Enterocyte",
                        "TA","Goblet" ,"Goblet" , "Goblet","Enteroendocrine",
                        "Paneth","Stem"  , "Stem",
                        "Stem","Stem","Tuft")

d10x.primary_exp_MHC<-d10x_exp_epi[,which(colnames(d10x_exp_epi)%in%MHCI)]
d10x.primary_exp_MHC$NAME<-rownames(d10x.primary_exp_MHC)
melt_exp<-melt(d10x.primary_exp_MHC)
plt_exp<-merge(plt, melt_exp, by="NAME")

plt_exp<-plt_exp[which(plt_exp$Site == "TI"),]


########### 
### count boxplot
###########

levels(plt_exp$Celltype)<-c("Entero-\nendocrine","BEST4\nEnterocyte","Enterocyte","TA","Goblet","Paneth","Stem","Tuft")
plt_exp$disease__ontology_label<-as.factor(plt_exp$disease__ontology_label)
levels(plt_exp$disease__ontology_label)<-c("CD","Control")

diff_exp_all_MHC<-diff_exp_sig_TI[which(diff_exp_sig_TI$gene%in%MHCI),]
diff_exp_all_MHC$celltype<-sapply(1:nrow(diff_exp_all_MHC), function(x) strsplit(diff_exp_all_MHC$cell.1[x],"_")[[1]][1])

highlight<-diff_exp_all_MHC[,c("gene","celltype","avg_log2FC")]
colnames(highlight)<-c("variable","Celltype","avg_log2FC")
highlight$direction<-sapply(1:nrow(highlight), function(x) if(highlight$avg_log2FC[x]<0){"down"}else{"up"})

highlight$Celltype<-as.factor(highlight$Celltype)
levels(highlight$Celltype)<-c("BEST4\nEnterocyte","Goblet","Paneth","Stem","TA","Tuft")

ggplot()+geom_violin(aes(disease__ontology_label, value, fill=disease__ontology_label),plt_exp, fill="grey", color="grey")+
  geom_boxplot(aes(disease__ontology_label, value, fill=disease__ontology_label),plt_exp, width=0.25, outlier.shape=NA)+
  facet_grid(variable~Celltype, scales="free_y")+theme_bw()+th_present+
  scale_fill_manual(values=c("dodgerblue3","lightgrey"))+
  theme(legend.position = "none")+ylab("Expression Level")+xlab("Diagnosis")+
  geom_rect(data = highlight, 
            fill = NA, aes(color=direction), xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, size=1.75)+
  scale_color_manual(values=c("green","yellow"))
ggsave(file=here("../EBI/scRNAseq_codon/figs","Kong_TI_MHC_singlecell.pdf"), w=8, h=10)
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_TI_MHC_singlecell.jpeg"), w=8, h=10)






########### 
### heatmap
###########
plt_epi<-merge(meta, TI_UMAP, by="NAME")

plt<-meta

plt<-plt[which(plt$Celltype != "Epithelial HBB HBA"),]
plt$Celltype<-as.factor(plt$Celltype)
levels(plt$Celltype)[which(levels(plt$Celltype)%in%unique(plt_epi$Celltype))]<-c("Enteroendocrine","BEST4 Enterocyte","Enterocyte","Enterocyte","Enterocyte",
                                                                                 "TA","Goblet" ,"Goblet" , "Goblet","Enteroendocrine",
                                                                                 "Paneth","Stem"  , "Stem",
                                                                                 "Stem","Stem","Tuft")

d10x.primary_exp_MHC<-d10x_exp[,which(colnames(d10x_exp)%in%MHCI)]
d10x.primary_exp_MHC$NAME<-rownames(d10x.primary_exp_MHC)
melt_exp<-melt(d10x.primary_exp_MHC)
plt_exp<-merge(plt, melt_exp, by="NAME")

plt_exp<-plt_exp[which(plt_exp$Site == "TI"),]

plt_exp$disease__ontology_label<-as.factor(plt_exp$disease__ontology_label)
levels(plt_exp$disease__ontology_label)<-c("CD","Control")


### scale each genen across cells then take mean for gene in each cell type
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_exp_scaled <- plt_exp %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(value))
plt_exp_summary <- plt_exp_scaled %>% 
  group_by(variable, Celltype, disease__ontology_label) %>%
  dplyr::summarize(Mean = mean(scaled, na.rm=TRUE))
plt_exp_summary<-as.data.frame(plt_exp_summary)


epi<-c("BEST4 Enterocyte","Stem","Enterocyte","Enteroendocrine",
       "Goblet", "Paneth","Tuft","TA")
plt_exp_summary$compartment<-"Non-epithelial"
plt_exp_summary$compartment[which(plt_exp_summary$Celltype%in%epi)]<-"Epithelial"

plt_exp_summary$label<-paste(plt_exp_summary$disease__ontology_label, plt_exp_summary$compartment, sep="\n")


## all cell types
ggplot(plt_exp_summary, aes(variable, Celltype, fill=Mean))+
  geom_tile()+facet_grid(label~., scales = "free_y", space = "free_y")+
  scale_fill_distiller(palette = "RdBu", name="Scaled\nMean\nExpression")+th+theme_classic()+
  ylab("")+xlab("")
ggsave(file="/home/redgar/Documents/EBI/scRNAseq_codon/figs/Kong_immunity_all_cells_MHCI_TI.pdf", w=9, h=15)




####
## select cell types
####
plt_exp_select<-plt_exp[which(plt_exp$Celltype %in% c("Mature DCs","Fibroblasts KCNN3 LY6H", epi)),]
plt_exp_select$Celltype<-factor(plt_exp_select$Celltype, levels=rev(c("Goblet", "Enterocyte","BEST4 Enterocyte","Enteroendocrine","Tuft",
                                                                      "TA","Paneth","Stem","Mature DCs","Fibroblasts KCNN3 LY6H")))

plt_exp_scaled_select <- plt_exp_select %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(value))
plt_exp_summary_select <- plt_exp_scaled_select %>% 
  group_by(variable, Celltype, disease__ontology_label) %>%
  dplyr::summarize(Mean = mean(scaled, na.rm=TRUE))
plt_exp_summary_select<-as.data.frame(plt_exp_summary_select)

epi<-c("BEST4 Enterocyte","Stem","Enterocyte","Enteroendocrine",
       "Goblet", "Paneth","Tuft","TA")
plt_exp_summary_select$compartment<-"Non-epithelial"
plt_exp_summary_select$compartment[which(plt_exp_summary_select$Celltype%in%epi)]<-"Epithelial"

plt_exp_summary_select$label<-paste(plt_exp_summary_select$disease__ontology_label, plt_exp_summary_select$compartment, sep="\n")


cell_count<-plt_exp[,c("NAME","disease__ontology_label","Celltype")]
cell_count<-cell_count[!duplicated(cell_count),]
cell_count$Celltype<-as.character(cell_count$Celltype)
total_count<-tapply(cell_count$NAME, list(cell_count$Celltype, cell_count$disease__ontology_label), function(x) length(unique(x)))
total_count<-melt(total_count)

total_count$compartment<-"Non-epithelial"
total_count$compartment[which(total_count$Var1%in%epi)]<-"Epithelial"

total_count$label<-paste(total_count$Var2, total_count$compartment, sep="\n")

## select cell types
total_count_select<-total_count[which(total_count$Var1 %in% c("Mature DCs","Fibroblasts KCNN3 LY6H", epi)),]
total_count_select$Celltype<-factor(total_count_select$Var1, levels=rev(c("Goblet", "Enterocyte","BEST4 Enterocyte","Enteroendocrine","Tuft",
                                                                          "TA","Paneth","Stem","Mature DCs","Fibroblasts KCNN3 LY6H")))


## add makr for significant
diff_exp_all_MHC<-diff_exp_sig_TI[which(diff_exp_sig_TI$gene%in%MHCI),]
diff_exp_all_MHC$celltype<-sapply(1:nrow(diff_exp_all_MHC), function(x) strsplit(diff_exp_all_MHC$cell.1[x],"_")[[1]][1])

highlight<-diff_exp_all_MHC[,c("gene","celltype","avg_log2FC")]
colnames(highlight)<-c("variable","Celltype","avg_log2FC")
highlight$direction<-sapply(1:nrow(highlight), function(x) if(highlight$avg_log2FC[x]<0){"down"}else{"up"})

highlight$Celltype<-as.factor(highlight$Celltype)

plt_exp_summary_select$sig<-sapply(1:nrow(plt_exp_summary_select), function(x){
  cell_sig<-diff_exp_all_MHC[which(diff_exp_all_MHC$celltype==plt_exp_summary_select$Celltype[x]),]
  if(nrow(cell_sig)==0){""}else{
    if(as.character(plt_exp_summary_select$variable[x])%in%cell_sig$gene){"*"}else{""}
  }
})

count_plt<-ggplot(total_count_select, aes(value, Celltype, fill=Var2, group=Var2))+
  geom_bar(stat = "identity", position=position_dodge(), color="black")+
  facet_grid(compartment~.,scales="free_y", space = "free_y")+th_present+theme_classic()+
  fillscale_diagnosis+xlab("Cell Count")+ylab("")+ theme(legend.position="top",
                                                         axis.title.y=element_blank(),
                                                         axis.text.y=element_blank(),
                                                         axis.ticks.y=element_blank(),
                                                         strip.background = element_blank(),
                                                         strip.text = element_blank(),
                                                         plot.margin = margin(1, 0.5, 0.25, 0, "cm"))

heat_plt<-ggplot(plt_exp_summary_select, aes(variable, Celltype, fill=Mean))+
  geom_tile(color="black")+facet_grid(compartment~disease__ontology_label, scales = "free_y", space = "free_y")+
  th_present+theme_classic()+
  ylab("")+xlab("")+ theme(legend.position="top",
                           strip.background.y = element_blank(),
                           strip.text.y = element_blank())+
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev = T, name ="Mean\nScaled\nExpression",
                                   p3 = .1, p4 = .6, p1 = .2, p2 = .6) +
  geom_text(aes(label=sig))

plot_grid(heat_plt, count_plt, rel_widths = c(1,0.175), align = "h", axis = "tb")

ggsave(file="/home/redgar/Documents/EBI/scRNAseq_codon/figs/Kong_immunity_TI_select_cells_MHCI_cellcount.pdf", w=16, h=5)





############
## Heat plot MHCI CO
############
load(here("../EBI/scRNAseq_codon/data/kong_immunity","Kong_Immunity_epithelial_MHCI_raw.RData"))
CO_UMAP<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/cluster/CO_EPI.scp.X_umap.coords.txt", header=T)
CO_UMAP<-CO_UMAP[-1,]
meta<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/metadata/scp_metadata_combined.v2.txt",header=T)


plt<-merge(meta, CO_UMAP, by="NAME")

plt$Celltype<-as.factor(plt$Celltype)
levels(plt$Celltype)<-c("BEST4 Enterocyte","Enterocyte",
                       "Enterocyte", "Enteroendocrine",
                       "TA","Goblet" ,"Goblet" , "Goblet",
                       "Paneth","Stem"  ,
                       "Stem","Stem","Tuft")

d10x.primary_exp_MHC<-d10x_exp[,which(colnames(d10x_exp)%in%MHCI)]
d10x.primary_exp_MHC$NAME<-rownames(d10x.primary_exp_MHC)
melt_exp<-melt(d10x.primary_exp_MHC)
plt_exp<-merge(plt, melt_exp, by="NAME")


plt_exp$disease__ontology_label<-as.factor(plt_exp$disease__ontology_label)
levels(plt_exp$disease__ontology_label)<-c("CD","Control")

########### 
### count boxplot
###########
plt_exp$Celltype<-factor(plt_exp$Celltype, levels=c("Enteroendocrine","BEST4 Enterocyte","Enterocyte","TA","Goblet","Paneth","Stem","Tuft"))
levels(plt_exp$Celltype)<-c("Entero-\nendocrine","BEST4\nEnterocyte","Enterocyte","TA","Goblet","Paneth","Stem","Tuft")
plt_exp$disease__ontology_label<-as.factor(plt_exp$disease__ontology_label)
levels(plt_exp$disease__ontology_label)<-c("CD","Control")

diff_exp_all_MHC<-diff_exp_sig_CO[which(diff_exp_sig_CO$gene%in%MHCI),]
diff_exp_all_MHC$celltype<-sapply(1:nrow(diff_exp_all_MHC), function(x) strsplit(diff_exp_all_MHC$cell.1[x],"_")[[1]][1])

highlight<-diff_exp_all_MHC[,c("gene","celltype","avg_log2FC")]
colnames(highlight)<-c("variable","Celltype","avg_log2FC")
highlight$direction<-sapply(1:nrow(highlight), function(x) if(highlight$avg_log2FC[x]<0){"down"}else{"up"})

highlight$Celltype<-as.factor(highlight$Celltype)
levels(highlight$Celltype)<-c("BEST4\nEnterocyte","Entero-\nendocrine","Goblet","Stem","TA")

ggplot()+geom_violin(aes(disease__ontology_label, value, fill=disease__ontology_label),plt_exp, fill="grey", color="grey")+
  geom_boxplot(aes(disease__ontology_label, value, fill=disease__ontology_label),plt_exp, width=0.25, outlier.shape=NA)+
  facet_grid(variable~Celltype, scales="free_y")+theme_bw()+th_present+
  scale_fill_manual(values=c("dodgerblue3","lightgrey"))+
  theme(legend.position = "none")+ylab("Expression Level")+xlab("Diagnosis")+
  geom_rect(data = highlight, 
            fill = NA, aes(color=direction), xmin = -Inf,xmax = Inf,
            ymin = -Inf,ymax = Inf, size=1.75)+
  scale_color_manual(values=c("green","yellow"))
ggsave(file=here("../EBI/scRNAseq_codon/figs","Kong_CO_MHC_singlecell.pdf"), w=8, h=10)
ggsave(file=here("../EBI/scRNAseq_codon/figs/jpeg","Kong_CO_MHC_singlecell.jpeg"), w=8, h=10)


################
## heat map
###############
plt_epi<-merge(meta, CO_UMAP, by="NAME")

plt<-meta

plt$Celltype<-as.factor(plt$Celltype)
levels(plt$Celltype)[which(levels(plt$Celltype)%in%unique(plt_epi$Celltype))]<-c("BEST4 Enterocyte","Enterocyte",
                                                                                 "Enterocyte", "Enteroendocrine",
                                                                                 "TA","Goblet" ,"Goblet" , "Goblet",
                                                                                 "Paneth","Stem"  ,
                                                                                 "Stem","Stem","Tuft")

d10x.primary_exp_MHC<-d10x_exp[,which(colnames(d10x_exp)%in%MHCI)]
d10x.primary_exp_MHC$NAME<-rownames(d10x.primary_exp_MHC)
melt_exp<-melt(d10x.primary_exp_MHC)
plt_exp<-merge(plt, melt_exp, by="NAME")

plt_exp<-plt_exp[which(plt_exp$Site == "TI"),]

plt_exp$disease__ontology_label<-as.factor(plt_exp$disease__ontology_label)
levels(plt_exp$disease__ontology_label)<-c("CD","Control")




### scale each genen across cells then take mean for gene in each cell type
scale_this <- function(x){(x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)}

plt_exp_scaled <- plt_exp %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(value))
plt_exp_summary <- plt_exp_scaled %>% 
  group_by(variable, Celltype, disease__ontology_label) %>%
  dplyr::summarize(Mean = mean(scaled, na.rm=TRUE))
plt_exp_summary<-as.data.frame(plt_exp_summary)


epi<-c("BEST4 Enterocyte","Stem","Enterocyte","Enteroendocrine",
       "Goblet", "Paneth","Tuft","TA")
plt_exp_summary$compartment<-"Non-epithelial"
plt_exp_summary$compartment[which(plt_exp_summary$Celltype%in%epi)]<-"Epithelial"

plt_exp_summary$label<-paste(plt_exp_summary$disease__ontology_label, plt_exp_summary$compartment, sep="\n")


## all cell types
ggplot(plt_exp_summary, aes(variable, Celltype, fill=Mean))+
  geom_tile()+facet_grid(label~., scales = "free_y", space = "free_y")+
  scale_fill_distiller(palette = "RdBu", name="Scaled\nMean\nExpression")+th+theme_classic()+
  ylab("")+xlab("")
ggsave(file="/home/redgar/Documents/EBI/scRNAseq_codon/figs/Kong_immunity_all_cells_MHCI_CO.pdf", w=9, h=15)




####
## select cell types
####
plt_exp_select<-plt_exp[which(plt_exp$Celltype %in% c("Mature DCs","Fibroblasts KCNN3 LY6H", epi)),]
plt_exp_select$Celltype<-factor(plt_exp_select$Celltype, levels=rev(c("Goblet", "Enterocyte","BEST4 Enterocyte","Enteroendocrine","Tuft",
                                                                      "TA","Paneth","Stem","Mature DCs","Fibroblasts KCNN3 LY6H")))

plt_exp_scaled_select <- plt_exp_select %>% group_by(variable) %>%
  dplyr::mutate(scaled = scale_this(value))
plt_exp_summary_select <- plt_exp_scaled_select %>% 
  group_by(variable, Celltype, disease__ontology_label) %>%
  dplyr::summarize(Mean = mean(scaled, na.rm=TRUE))
plt_exp_summary_select<-as.data.frame(plt_exp_summary_select)

epi<-c("BEST4 Enterocyte","Stem","Enterocyte","Enteroendocrine",
       "Goblet", "Paneth","Tuft","TA")
plt_exp_summary_select$compartment<-"Non-epithelial"
plt_exp_summary_select$compartment[which(plt_exp_summary_select$Celltype%in%epi)]<-"Epithelial"

plt_exp_summary_select$label<-paste(plt_exp_summary_select$disease__ontology_label, plt_exp_summary_select$compartment, sep="\n")


cell_count<-plt_exp[,c("NAME","disease__ontology_label","Celltype")]
cell_count<-cell_count[!duplicated(cell_count),]
cell_count$Celltype<-as.character(cell_count$Celltype)
total_count<-tapply(cell_count$NAME, list(cell_count$Celltype, cell_count$disease__ontology_label), function(x) length(unique(x)))
total_count<-melt(total_count)

total_count$compartment<-"Non-epithelial"
total_count$compartment[which(total_count$Var1%in%epi)]<-"Epithelial"

total_count$label<-paste(total_count$Var2, total_count$compartment, sep="\n")

## select cell types
total_count_select<-total_count[which(total_count$Var1 %in% c("Mature DCs","Fibroblasts KCNN3 LY6H", epi)),]
total_count_select$Celltype<-factor(total_count_select$Var1, levels=rev(c("Goblet", "Enterocyte","BEST4 Enterocyte","Enteroendocrine","Tuft",
                                                                          "TA","Paneth","Stem","Mature DCs","Fibroblasts KCNN3 LY6H")))


## add makr for significant
diff_exp_all_MHC<-diff_exp_sig_CO[which(diff_exp_sig_CO$gene%in%MHCI),]
diff_exp_all_MHC$celltype<-sapply(1:nrow(diff_exp_all_MHC), function(x) strsplit(diff_exp_all_MHC$cell.1[x],"_")[[1]][1])


plt_exp_summary_select$sig<-sapply(1:nrow(plt_exp_summary_select), function(x){
  cell_sig<-diff_exp_all_MHC[which(diff_exp_all_MHC$celltype==plt_exp_summary_select$Celltype[x]),]
  if(nrow(cell_sig)==0){""}else{
    if(as.character(plt_exp_summary_select$variable[x])%in%cell_sig$gene){"*"}else{""}
  }
})

count_plt<-ggplot(total_count_select, aes(value, Celltype, fill=Var2, group=Var2))+
  geom_bar(stat = "identity", position=position_dodge(), color="black")+
  facet_grid(compartment~.,scales="free_y", space = "free_y")+th_present+theme_classic()+
  fillscale_diagnosis+xlab("Cell Count")+ylab("")+ theme(legend.position="top",
                                                         axis.title.y=element_blank(),
                                                         axis.text.y=element_blank(),
                                                         axis.ticks.y=element_blank(),
                                                         strip.background = element_blank(),
                                                         strip.text = element_blank(),
                                                         plot.margin = margin(1, 0.5, 0.25, 0, "cm"))

heat_plt<-ggplot(plt_exp_summary_select, aes(variable, Celltype, fill=Mean))+
  geom_tile(color="black")+facet_grid(compartment~disease__ontology_label, scales = "free_y", space = "free_y")+
  th_present+theme_classic()+
  ylab("")+xlab("")+ theme(legend.position="top",
                           strip.background.y = element_blank(),
                           strip.text.y = element_blank())+
  scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0, rev = T, name ="Mean\nScaled\nExpression",
                                   p3 = .1, p4 = .6, p1 = .2, p2 = .6) +
  geom_text(aes(label=sig))

plot_grid(heat_plt, count_plt, rel_widths = c(1,0.175), align = "h", axis = "tb")

ggsave(file="/home/redgar/Documents/EBI/scRNAseq_codon/figs/Kong_immunity_CO_select_cells_MHCI_cellcount.pdf", w=16, h=5)






############
## T cell score correlation
###########
load(here("../EBI/scRNAseq_codon/data/kong_immunity","Kong_Immunity_epithelial_MHCI_raw.RData"))

meta<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/metadata/scp_metadata_combined.v2.txt",header=T)
meta <- meta[-(1),]

cell_type_count<-as.data.frame(table(meta$donor_id,meta$Site, meta$Celltype))



meta<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/metadata/scp_metadata_combined.v2.txt",header=T)
TI_UMAP<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/cluster/TI_EPI.scp.X_umap.coords.txt", header=T)
TI_UMAP<-TI_UMAP[-1,]
CO_UMAP<-read.table("/media/redgar/Seagate Portable Drive/kong_immunity/SCP1884/cluster/CO_EPI.scp.X_umap.coords.txt", header=T)
CO_UMAP<-CO_UMAP[-1,]

plt_TI<-merge(meta, TI_UMAP, by="NAME")
plt_TI<-merge(plt_TI, score_data_allcompartment, by="NAME")
plt_TI<-plt_TI[which(plt$Site == "TI"),]
plt_TI<-plt_TI[which(plt_TI$Celltype != "Epithelial HBB HBA"),]

mn_sample_MHCI_TI<-as.data.frame(tapply(plt_TI$MHCI_score1_allcompartment, plt_TI$donor_id, mean))
colnames(mn_sample_MHCI_TI)<-"Mean_MHCI_score"
mn_sample_MHCI_TI$donor_id<-rownames(mn_sample_MHCI_TI)
mn_sample_MHCI_TI$Site<-"TI"

plt_CO<-merge(meta, CO_UMAP, by="NAME")
plt_CO<-merge(plt_CO, score_data_allcompartment, by="NAME")

mn_sample_MHCI_CO<-as.data.frame(tapply(plt_CO$MHCI_score1_allcompartment, plt_CO$donor_id, mean))
colnames(mn_sample_MHCI_CO)<-"Mean_MHCI_score"
mn_sample_MHCI_CO$donor_id<-rownames(mn_sample_MHCI_CO)
mn_sample_MHCI_CO$Site<-"CO"

mn_sample_MHCI<-rbind(mn_sample_MHCI_CO, mn_sample_MHCI_TI)


cell_count_score<-merge(cell_type_count, mn_sample_MHCI,  by.x=c("Var1","Var2"),by.y=c("donor_id","Site"),)

meta_min<-meta[,c("donor_id","disease__ontology_label")]
meta_min<-meta_min[!duplicated(meta_min),]

cell_count_score<-merge(cell_count_score, meta_min, by.x="Var1", by.y="donor_id")



summary_stat<-as.data.frame(cell_count_score[which(cell_count_score$Var3%in%c("Tregs","T cells CD4 IL17A", "T cells CD4 FOSB" , "T cells CD8"   ,
                                                                              "T cells OGT","T cells CD8 KLRG1", "T cells Naive CD4")),] %>%  
                              group_by(Var3,Var2, disease__ontology_label) %>% summarise(cor_coef = cor.test(Mean_MHCI_score, Freq, method="spearman")$estimate,
                                                                       p_val = cor.test(Mean_MHCI_score, Freq, method="spearman")$p.value))


ggplot(cell_count_score[which(cell_count_score$Var3%in%c("Tregs","T cells CD4 IL17A", "T cells CD4 FOSB" , "T cells CD8"   ,
                                                         "T cells OGT","T cells CD8 KLRG1", "T cells Naive CD4")),], 
       aes(Mean_MHCI_score, Freq, fill=disease__ontology_label, color=disease__ontology_label))+
  scale_fill_manual(values=c("dodgerblue3","lightgrey"), name="Diagnosis")+
  scale_color_manual(values=c("dodgerblue3","lightgrey"), name="Diagnosis")+
  geom_point(shape=21, color="black")+facet_grid(Var3~Var2, scales="free_y")+theme_bw()+th_present+theme(legend.position = "none")+
  stat_smooth(method="lm",se=F)+
  geom_text(aes(x=0.5, y=500, label=paste("Rs = ",signif(cor_coef,2),";","p = ",signif(p_val,1))),data=summary_stat[which(summary_stat$disease__ontology_label=="normal"),])+
  geom_text(aes(x=0.5, y=800, label=paste("Rs = ",signif(cor_coef,2),";","p = ",signif(p_val,1))),data=summary_stat[which(summary_stat$disease__ontology_label=="Crohn's disease"),])
  


summary_stat<-as.data.frame(cell_count_score[which(cell_count_score$Var2%in%c("Activated T","CD4 T cell","CD8 T cell")),] %>%  
                              group_by(Var2) %>% summarise(cor_coef = cor.test(Mean_MHCI_score, Freq, method="spearman")$estimate,
                                                           p_val = cor.test(Mean_MHCI_score, Freq, method="spearman")$p.value))

ggplot(cell_count_score[which(cell_count_score$Var2%in%c("Activated T","CD4 T cell","CD8 T cell")),], 
       aes(Mean_MHCI_score, Freq))+
  geom_point(aes( fill=orig.ident),color="black",shape=21)+facet_wrap(~Var2, scales="free_y")+
  stat_smooth(method="lm",se=F, color="black")+
  fillscale_diagnosis+colscale_diagnosis+theme_bw()+th_present+theme(legend.position = "none")+
  geom_text(aes(x=0.3, y=10, label=paste("Rs = ",signif(cor_coef,2),";","p = ",signif(p_val,1))),data=summary_stat, color="black")


