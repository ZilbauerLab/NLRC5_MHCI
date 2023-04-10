
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
crypt_villis = c("CEACAM7", "PLAC8", "CEACAM1", "TSPAN1", "CEACAM5", "CEACAM6", "IFI27", "DHRS9", "KRT20", "RHOC", "CD177", "PKIB", "HPGD", "LYPD8", "APOBEC1", "APOB", "APOA4", "APOA1", "NPC1L1", "EGFR", "KLF4", "ENPP3", "NT5E", "SLC28A2", "ADA")
#"SEPP1" is not in the dataset


## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))

# no neonatal or posneg
d10x.primary<-subset(d10x.primary, subset = individual %in% c(
  "T017","T019","T176","T189","T197","T202","T203","T024","T036","T44","T057",
  "T160","T161","T175","T182","T184","T180"))

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
# no neonatal or posneg
primary.epi.cells<-subset(primary.epi.cells, subset = individual %in% c(
  "T017","T019","T176","T189","T197","T202","T203","T024","T036","T44","T057",
  "T160","T161","T175","T182","T184","T180"))


#filter remaining immune (179 cells)
primary.epi.cells<-subset(primary.epi.cells, subset = cluster_ID != "B cell" & cluster_ID != "T cell")


score_data_epi<-score_data[which(rownames(score_data)%in%colnames(primary.epi.cells)),]
full_geneset_correlation<-suppressWarnings(cor.test(score_data_epi$MHCI_score1,score_data_epi$crypt_villis_score1, method="spearman"))$estimate

print("Correaltion to report")
print(full_geneset_correlation)


score_jackknife<-lapply(1:length(c(MHCI,crypt_villis)),function(x){
  if(x<=length(MHCI)){
    gene<-MHCI[x]
    MHCI_jack<-MHCI[-x]
    crypt_villis_jack<-crypt_villis
  }else{
    gene<-crypt_villis[x-length(MHCI)]
    MHCI_jack<-MHCI
    crypt_villis_jack<-crypt_villis[-(x-length(MHCI))]
      }

  d10x.primary <- AddModuleScore(
    object = d10x.primary,
    features = list(MHCI_jack),
    ctrl = 5,
    name = 'MHCI_score'
  )
  
  d10x.primary <- AddModuleScore(
    object = d10x.primary,
    features = list(crypt_villis_jack),
    ctrl = 5,
    name = 'crypt_villis_score'
  )
  
  score_data<-d10x.primary@meta.data[,c("MHCI_score1","crypt_villis_score1")]
  score_data_epi<-score_data[which(rownames(score_data)%in%colnames(primary.epi.cells)),]
  
  data.frame(gene=gene, correlation=suppressWarnings(cor.test(score_data_epi$MHCI_score1,score_data_epi$crypt_villis_score1, method="spearman"))$estimate)})

score_jackknife<-do.call(rbind,score_jackknife)

score_jackknife$list<-"MHCI Genes"
score_jackknife$list[which(score_jackknife$gene%in%crypt_villis)]<-"Crypt Villus Genes"

ggplot(score_jackknife, aes(correlation))+geom_histogram(bins=10, color="black",fill="lightgrey")+
  geom_vline(xintercept = full_geneset_correlation)+theme_bw()+th+
  facet_wrap(~list)+xlab("Correlation")

ggsave(file=here("figs/jpeg","jack_knife_score_correlation.jpeg"), w=6, h=3)
ggsave(file=here("figs","jack_knife_score_correlation.pdf"), w=6, h=3)


score_jackknife$deviation<-full_geneset_correlation-score_jackknife$correlation
score_jackknife$deviation_clear<-round(score_jackknife$deviation,3)

score_jackknife[order(score_jackknife$deviation),]
