
#'---
#'title: MHCI and Crypt-villus score correlation
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
d10x.stim <- AddMetaData(d10x.stim, metadata = cell_typelabel)

## LogNormalize: Feature counts for each cell are divided by the total counts for that cell and multiplied by the scale.factor. This is then natural-log transformed using log1p.
# This is log(TP10K+1)
d10x.stim <- NormalizeData(d10x.stim,scale.factor = 10000, normalization.method = "LogNormalize")
d10x.stim_scores<-FetchData(object = d10x.stim, vars = c(MHCI,crypt_villis))


## diagional
cormat<-cor(d10x.stim_scores, method="spearman")

fill_mat<-cormat
fill_mat[upper.tri(fill_mat)]<- NA
fill_melted_cormat <- melt(fill_mat, na.rm = TRUE)

## remove same gene cor
fill_melted_cormat<-fill_melted_cormat[which(!(fill_melted_cormat$Var1==fill_melted_cormat$Var2)),]

## select some correlations to label
label_max_offset<-fill_melted_cormat[which(fill_melted_cormat$Var1%in%crypt_villis & fill_melted_cormat$Var2%in%MHCI),]
label_max_offset<-label_max_offset[!duplicated(label_max_offset),]

label_NLRC5<-fill_melted_cormat[which(fill_melted_cormat$Var1=="NLRC5" & fill_melted_cormat$Var2%in%MHCI),]

selected_labels<-rbind(label_NLRC5,label_max_offset[rev(order(label_max_offset$value)),][1:10,])

fill_melted_cormat$Var1<-factor(fill_melted_cormat$Var1, levels=c(MHCI,crypt_villis))
fill_melted_cormat$Var2<-factor(fill_melted_cormat$Var2, levels=c(MHCI,crypt_villis))


ggplot()+geom_tile(aes(Var1, Var2, fill=value),fill_melted_cormat)+
  geom_text(aes(Var1, Var2, label=round(value,2)), selected_labels, color="grey40", size=2.5)+
  scale_fill_gradient2(high="#2171b5",mid="#f3f7fc",low="red",
                       na.value="grey90", midpoint=0)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size =9, color="black"),
        legend.text = element_text(size =10),
        legend.title = element_text(size =11),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_y_discrete(position = "right") +
  geom_rect(aes(xmin=length(MHCI)-0.5, xmax=length(c(MHCI,crypt_villis))-1.5,ymax=length(MHCI)+0.5, ymin=0.5),fill=NA, color="grey")

ggsave(file=here("figs/jpeg","score_correlation.jpeg"), w=10, h=9)
ggsave(file=here("figs","score_correlation.pdf"), w=10, h=9)

 

#################
#'## Are these correlations more than expected by chance?
#################
label_max_offset<-fill_melted_cormat[which(fill_melted_cormat$Var1%in%crypt_villis & fill_melted_cormat$Var2%in%MHCI),]
median(abs(label_max_offset$value))



#'###compare MHCI to random
suppressWarnings(random_medians_MHCtornd<-sapply(1:10000, function(x){
  set.seed(x)
  ## random gene list
  rnd<-sample(rownames(d10x.stim),length(crypt_villis))
  
  ## fetch random gene list
  d10x.stim_MHCI_rnd<-FetchData(object = d10x.stim, vars = c(MHCI,rnd))
  
  ## diagional
  cormat<-cor(d10x.stim_MHCI_rnd, method="spearman")
  
  fill_mat_rnd<-cormat
  fill_mat_rnd[upper.tri(fill_mat_rnd)]<- NA
  fill_melted_cormat_rnd <- melt(fill_mat_rnd, na.rm = TRUE)
  
  ## remove same gene cor
  fill_melted_cormat_rnd<-fill_melted_cormat_rnd[which(!(fill_melted_cormat_rnd$Var1==fill_melted_cormat_rnd$Var2)),]
  
  ## select some correlations to label
  label_max_offset_rnd<-fill_melted_cormat_rnd[which(fill_melted_cormat_rnd$Var1%in%rnd & fill_melted_cormat_rnd$Var2%in%MHCI),]
  median(abs(label_max_offset_rnd$value))
}))


## how many larger
length(which(random_medians_MHCtornd>=median(abs(label_max_offset$value))))
ggplot()+geom_histogram(aes(random_medians_MHCtornd))+geom_vline(aes(xintercept=median(abs(label_max_offset$value))))

(length(which(random_medians_MHCtornd>=median(abs(label_max_offset$value))))+1)/(10000+1)






#'### compare crypt villus to random
suppressWarnings(random_medians_crypviltornd<-sapply(1:10000, function(x){
  set.seed(x)
  ## random gene list
  rnd<-sample(rownames(d10x.stim),length(MHCI))
  
  ## fetch random gene list
  d10x.stim_crypt_villis_rnd<-FetchData(object = d10x.stim, vars = c(rnd,crypt_villis))
  
  ## diagional
  cormat<-cor(d10x.stim_crypt_villis_rnd, method="spearman")
  
  fill_mat_rnd<-cormat
  fill_mat_rnd[upper.tri(fill_mat_rnd)]<- NA
  fill_melted_cormat_rnd <- melt(fill_mat_rnd, na.rm = TRUE)
  
  ## remove same gene cor
  fill_melted_cormat_rnd<-fill_melted_cormat_rnd[which(!(fill_melted_cormat_rnd$Var1==fill_melted_cormat_rnd$Var2)),]
  
  ## select some correlations to label
  label_max_offset_rnd<-fill_melted_cormat_rnd[which(fill_melted_cormat_rnd$Var1%in%crypt_villis & fill_melted_cormat_rnd$Var2%in%rnd),]
  median(abs(label_max_offset_rnd$value))
}))


## how many larger
length(which(random_medians_crypviltornd>=median(abs(label_max_offset$value))))
ggplot()+geom_histogram(aes(random_medians_crypviltornd))+geom_vline(aes(xintercept=median(abs(label_max_offset$value))))

(length(which(random_medians_crypviltornd>=median(abs(label_max_offset$value))))+1)/(10000+1)




                  
#######################                  
#'## Individual correlation plots
#######################


gene_correlation<-function(gene_name1, gene_name2){
  plt_gene<-FetchData(object = d10x.stim, vars = c(gene_name1, gene_name2))
  #plt_gene$cell<-rownames(plt_gene)

  ## ignore 0s for correlation?
  plt_gene$zero<-"Both non-zero"
  plt_gene$zero[which(plt_gene[,1]==0 | plt_gene[,2]==0)]<-"One or both\ngenes zero"
  
  allcor<-cor(plt_gene[,1],plt_gene[,2], method="spearman")
  #non-zero cor
  plt_gene_non_zero<-plt_gene[which(plt_gene$zero=="Both non-zero"),]
  nonzerocor<-cor(plt_gene_non_zero[,1],plt_gene_non_zero[,2], method="spearman")
  
  scatter<-ggplot()+geom_point(aes(plt_gene[,1], plt_gene[,2],  color=zero),plt_gene)+ylab(gene_name2)+xlab(gene_name1)+
    theme_bw()+th+scale_color_manual(values=c("cornflowerblue","lightgrey"), name="Both genes\ndetected in cell")+
    stat_smooth(aes(plt_gene[,1], plt_gene[,2]),plt_gene ,method="lm", se=F, color="black")+#
    stat_smooth(aes(plt_gene_non_zero[,1], plt_gene_non_zero[,2]),plt_gene_non_zero ,method="lm", se=F, color="#1f66e5")+
    theme(legend.text=element_text(size=8),
          legend.title=element_text(size=10))
  
  get_leg = function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend}
  
  leg_scatter = get_leg(scatter)
  
  
  density1<-ggplot(plt_gene, aes(plt_gene[,1]))+geom_density(fill="#a6bddb")+
    theme_void()+th+xlab(gene_name1)+theme(axis.title=element_blank(),
                                         axis.text=element_blank(),
                                         axis.ticks=element_blank(),
                                         legend.position = "none",
                                         panel.background = element_blank(),
                                         panel.grid.major = element_blank(), 
                                         panel.grid.minor = element_blank(),
                                         plot.margin = unit(c(0.75,0.2,-0.1,1), "cm"))
  
  density2<-ggplot(plt_gene, aes(plt_gene[,2]))+geom_density(fill="#a6bddb")+
    coord_flip()+theme_void()+th+xlab(gene_name1)+theme(axis.title=element_blank(),
                                                      axis.text=element_blank(),
                                                      axis.ticks=element_blank(),
                                                      legend.position = "none",
                                                      panel.background = element_blank(),
                                                      panel.grid.major = element_blank(), 
                                                      panel.grid.minor = element_blank(),
                                                      plot.margin = unit(c(0.2,1,1.1,-0.1), "cm"))
  
  placeholder<-ggplot() + theme_void()
  
  corvalues<-ggplot(plt_gene,aes(plt_gene[,1], plt_gene[,2],  color=zero)) + theme_void() +
    annotate("text", x=max(plt_gene[,1])*0.25, y=max(plt_gene[,2])*0.5, label = paste("All cells\nR =",round(allcor,2)))+
    annotate("text", x=max(plt_gene[,1])*0.75, y=max(plt_gene[,2])*0.5, label = paste("Both genes non-zero\nR =",round(nonzerocor,2)), color="#1f66e5")+
    xlim(0,max(plt_gene[,1]))+ylim(0,max(plt_gene[,2]))
  
  grid_fig<-plot_grid(density1, leg_scatter, 
                      scatter+theme(legend.position = "none"), density2, 
                      corvalues,placeholder,
                      ncol = 2, axis="lr", 
                      rel_heights = c(1,4,0.5),
                      rel_widths = c(4,1))
  
  ggsave(file=paste(here("figs/jpeg/"),gene_name1,"_",gene_name2,"_correlation.jpeg", sep=""), w=7, h=7)
  ggsave(file=paste(here("figs/"),gene_name1,"_",gene_name2,"_correlation.pdf", sep=""), w=7, h=7)
  
  grid_fig
}

gene_correlation("NLRC5","TAP1")
gene_correlation("NLRC5","PSMB9")
gene_correlation("B2M","HLA-B")

gene_correlation("PLAC8","B2M")
gene_correlation("IFI27","HLA-A")



######### 
#'## summary statistics on correlation
######### 
gene_combos<-combn(MHCI, 2)

gene_correlation_values<-do.call(rbind, lapply(1:ncol(gene_combos), function(x){
  gene_name1=gene_combos[1,x]
  gene_name2=gene_combos[2,x]
  plt_gene<-FetchData(object = d10x.stim, vars = c(gene_name1, gene_name2))
  
  ## ignore 0s for correlation
  plt_gene$zero<-"Both non-zero"
  plt_gene$zero[which(plt_gene[,1]==0 | plt_gene[,2]==0)]<-"One or both\ngenes zero"
  
  allcor<-cor(plt_gene[,1],plt_gene[,2], method="spearman")
  #non-zero cor
  plt_gene_non_zero<-plt_gene[which(plt_gene$zero=="Both non-zero"),]
  nonzerocor<-cor(plt_gene_non_zero[,1],plt_gene_non_zero[,2], method="spearman")
  data.frame(gene_name1=gene_name1, gene_name2=gene_name2, nonzerocor=nonzerocor, allcor=allcor)
}))
write.csv(gene_correlation_values, paste(here("output/"),"organoid_MHC1_correlation_edges.csv", sep=""), row.names=F, quote = F)



plt_gene_MHCI<-FetchData(object = d10x.stim, vars = MHCI)
counts_cells<-apply(plt_gene_MHCI, 2, function(x) length(which(x!=0)))
mean_cells<-apply(plt_gene_MHCI, 2, mean)

node_values<-cbind(melt(counts_cells), melt(mean_cells))
colnames(node_values)<-c("cell_expressing","mean_expression")
node_values$gene=rownames(node_values)

write.csv(node_values, paste(here("output/"),"organoid_MHC1_correlation_node.csv", sep=""), row.names=F, quote = F)


##### MHCI correlation only in non zero alternate to network
gene_correlation_values$gene_name1<-factor(gene_correlation_values$gene_name1, levels=rev(MHCI))
gene_correlation_values$gene_name2<-factor(gene_correlation_values$gene_name2, levels=rev(MHCI))

selected_labels<-gene_correlation_values[which(gene_correlation_values$gene_name2=="NLRC5"),]

ggplot()+geom_tile(aes(gene_name1, gene_name2, fill=nonzerocor),gene_correlation_values, color="black")+
  geom_text(aes(gene_name1, gene_name2, label=round(nonzerocor,2)), selected_labels, color="grey10", size=4)+
  scale_fill_gradient2(high="#2171b5",mid="#f3f7fc",low="red",
                       na.value="grey90", midpoint=0.1, name="Correlation", breaks=seq(-1,1,0.2), limits=c(-1,1))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_text(size =12, color="black"),
        legend.text = element_text(size =10),
        legend.title = element_text(size =11),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  scale_y_discrete(position = "right") 

ggsave(file=here("figs/jpeg","MHCI_score_correlation_organoid_nonzero.jpeg"), w=10, h=9)
ggsave(file=here("figs","MHCI_score_correlation_organoid_nonzero.pdf"), w=10, h=9)



#'## R Session Info
sessionInfo()
