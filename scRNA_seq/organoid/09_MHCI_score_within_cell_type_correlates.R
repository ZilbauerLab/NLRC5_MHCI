

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
cell_typelabel$cluster_ID[which(cell_typelabel$cluster_ID=="entero_stem")]<-"enterocyte_stem"
d10x.stim <- AddMetaData(d10x.stim, metadata = cell_typelabel)

##LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.stim <- NormalizeData(d10x.stim,scale.factor = 10000, normalization.method = "LogNormalize")

hsp<-rownames(d10x.stim)[grep("^HSP",rownames(d10x.stim))]

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

d10x.stim <- AddModuleScore(
  object = d10x.stim,
  features = list(hsp),
  ctrl = 5,
  name = 'HSP_score'
)


d10x.stim@meta.data$cluster_ID<-factor(d10x.stim@meta.data$cluster_ID, levels=c("crypt","TA","enterocyte_mt","enterocyte"))

#score_data_QC<-d10x.stim@meta.data[,c("MHCI_score1","crypt_villis_score1","cluster_ID","")]
# ggplot(d10x.stim@meta.data, aes(percent.mt, MHCI_score1))+geom_point()+stat_smooth(method="lm", se=F)+
#   facet_grid(orig.ident~cluster_ID)

#' within a cell type
grouped<-d10x.stim@meta.data %>% group_by(cluster_ID, orig.ident) 
correlations<-as.data.frame(summarize(grouped, cor(MHCI_score1,percent.mt, method="spearman")))
colnames(correlations)[3]<-"correlation"

ggplot(d10x.stim@meta.data, aes(MHCI_score1, percent.mt))+
  geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
  facet_grid(orig.ident~cluster_ID)+scale_fill_manual(values=c("#6baed6","#f16913","#78c679","#238b45"))+
  geom_label(data=correlations, aes(x=Inf, y=Inf, label = paste("R = ", round(correlation,2), sep = " ")), hjust = 1, vjust = 1)
ggsave(file=here("figs/jpeg/MHCI_percentmt.jpeg"), w=10, h=6)
ggsave(file=here("figs/MHCI_percentmt.pdf"), w=10, h=6)



#' across all cells (could be driven my entero mt)
ggplot(d10x.stim@meta.data, aes(percent.mt, MHCI_score1))+geom_point()+stat_smooth(method="lm", se=F)+
  facet_grid(.~orig.ident)

#' across all cells (excluding entero mt)
ggplot(d10x.stim@meta.data[which(d10x.stim@meta.data$cluster_ID!="enterocyte_mt"),], aes(percent.mt, MHCI_score1))+
  geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
  facet_grid(.~orig.ident)+scale_fill_manual(values=c("#6baed6","#f16913","#238b45"))

#' crypt villus also higher in high mt?
ggplot(d10x.stim@meta.data[which(d10x.stim@meta.data$cluster_ID!="enterocyte_mt"),], aes(percent.mt, crypt_villis_score1))+geom_point()+stat_smooth(method="lm", se=F)+
facet_grid(orig.ident~cluster_ID)

ggplot(d10x.stim@meta.data[which(d10x.stim@meta.data$cluster_ID!="enterocyte_mt"),], aes(percent.mt, crypt_villis_score1))+
  geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
  facet_grid(orig.ident~.)+  scale_fill_manual(values=c("#6baed6","#f16913","#238b45"))


#' # Percent mt should be higher in IFNg over all then
#' Well you sort of alrady say this cause you filter more cells from both treatments then controls


#' Are the crypt with e higher MHCI more differentiated? No
grouped<-d10x.stim@meta.data %>% group_by(cluster_ID, orig.ident) 
correlations<-as.data.frame(summarize(grouped, cor(MHCI_score1,crypt_villis_score1, method="spearman")))
colnames(correlations)[3]<-"correlation"

ggplot(d10x.stim@meta.data, aes(MHCI_score1, crypt_villis_score1))+
  geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
  facet_grid(orig.ident~cluster_ID)+scale_fill_manual(values=c("#6baed6","#f16913","#238b45"))+
  geom_label(data=correlations, aes(x=Inf, y=Inf, label = paste("R = ", round(correlation,2), sep = " ")), hjust = 1, vjust = 1)

ggsave(file=here("figs/jpeg/MHCI_crypt_villus.jpeg"), w=10, h=6)
ggsave(file=here("figs/MHCI_crypt_villus.pdf"), w=10, h=6)



#' Are the cells more stressed out (HSP score)
ggplot(d10x.stim@meta.data, aes(MHCI_score1, HSP_score1))+
  geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
  facet_grid(orig.ident~cluster_ID)+scale_fill_manual(values=c("#6baed6","#f16913","#238b45"))




##################
#'## Stress genes with increasing MHCI score in IFNg enterocytes
##################
d10x.stim_scores<-FetchData(object = d10x.stim, vars = c("DNAJB1","HSPA1A"))
d10x.stim_scores$ID<-rownames(d10x.stim_scores)
d10x.stim@meta.data$ID<-rownames(d10x.stim@meta.data)

plt<-merge(d10x.stim@meta.data, d10x.stim_scores, by="ID")

plt$DNAJB1_zero<-"zero"
plt$DNAJB1_zero[which(plt$DNAJB1!=0)]<-"non zero"

plt$HSPA1A_zero<-"zero"
plt$HSPA1A_zero[which(plt$HSPA1A!=0)]<-"non zero"


# ggplot(plt, aes(MHCI_score1, DNAJB1))+
#   geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
#   facet_grid(cluster_ID~orig.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))


grouped<-plt %>% group_by(cluster_ID, orig.ident, DNAJB1_zero) 
correlations<-as.data.frame(summarize(grouped, cor(MHCI_score1,DNAJB1, method="spearman")))
colnames(correlations)[4]<-"correlation"
# remove correlations for not measured gene
correlations<-correlations[which(!(is.na(correlations$correlation))),]

ggplot(plt, aes(MHCI_score1, DNAJB1, group=DNAJB1_zero))+
  geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
  facet_grid(orig.ident~cluster_ID)+scale_fill_manual(values=c("#6baed6","#f16913","#78c679","#238b45"))+
  geom_label(data=correlations, aes(x=Inf, y=Inf, label = paste("R = ", round(correlation,2), sep = " ")), hjust = 1, vjust = 1)

ggsave(file=here("figs/jpeg/MHCI_DNAJB1.jpeg"), w=10, h=6)
ggsave(file=here("figs/MHCI_DNAJB1.pdf"), w=10, h=6)


grouped<-plt %>% group_by(cluster_ID, orig.ident, HSPA1A_zero) 
correlations<-as.data.frame(summarize(grouped, cor(MHCI_score1,HSPA1A, method="spearman")))
colnames(correlations)[4]<-"correlation"
# remove correlations for not measured gene
correlations<-correlations[which(!(is.na(correlations$correlation))),]

ggplot(plt, aes(MHCI_score1, HSPA1A, group=HSPA1A_zero))+
  geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
  facet_grid(orig.ident~cluster_ID)+scale_fill_manual(values=c("#6baed6","#f16913","#78c679","#238b45"))+
  geom_label(data=correlations, aes(x=Inf, y=Inf, label = paste("R = ", round(correlation,2), sep = " ")), hjust = 1, vjust = 1)

ggsave(file=here("figs/jpeg/MHCI_HSPA1A.jpeg"), w=10, h=6)
ggsave(file=here("figs/MHCI_HSPA1A.pdf"), w=10, h=6)

###############
#### IFNg receptors
###############
d10x.stim_scores<-FetchData(object = d10x.stim, vars = c("IFNGR1","IFNGR2"))
d10x.stim_scores$ID<-rownames(d10x.stim_scores)
d10x.stim@meta.data$ID<-rownames(d10x.stim@meta.data)

plt<-merge(d10x.stim@meta.data, d10x.stim_scores, by="ID")

plt$IFNGR1_zero<-"zero"
plt$IFNGR1_zero[which(plt$IFNGR1!=0)]<-"non zero"

plt$IFNGR2_zero<-"zero"
plt$IFNGR2_zero[which(plt$IFNGR2!=0)]<-"non zero"


ggplot(plt, aes(IFNGR1,MHCI_score1))+
  geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
  facet_grid(cluster_ID~orig.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))

ggplot(plt, aes(IFNGR1,MHCI_score1, group=IFNGR1_zero))+
  geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
  facet_grid(cluster_ID~orig.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))


ggplot(plt, aes( IFNGR2,MHCI_score1))+
  geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
  facet_grid(cluster_ID~orig.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))

ggplot(plt, aes( IFNGR2, MHCI_score1, group=IFNGR2_zero))+
  geom_point(aes(fill=cluster_ID),shape=21, color="black")+stat_smooth(method="lm", se=F)+
  facet_grid(cluster_ID~orig.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))



##################
#'## cell cycle?
##################
d10x.stim_norm<-readRDS(here("data","d10x_normalized.rds"))


d10x.stim_norm@meta.data$ID<-rownames(d10x.stim_norm@meta.data)

plt<-merge(d10x.stim@meta.data, d10x.stim_norm@meta.data[,c("ID","S.Score","G2M.Score","Phase")], by="ID")


ggplot(plt, aes(Phase,MHCI_score1))+
  geom_boxplot(aes(fill=cluster_ID),color="black")+stat_smooth(method="lm", se=F)+theme_bw()+th+
  facet_grid(cluster_ID~orig.ident)+scale_fill_manual(values=c("#6baed6","#238b45","#78c679","#f16913"))

ggplot(plt, aes(MHCI_score1))+
  geom_density(aes(fill=Phase),alpha=0.25,color="black")+theme_bw()+th+
  facet_grid(cluster_ID~orig.ident)
