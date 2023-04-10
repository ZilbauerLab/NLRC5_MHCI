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
library(cowplot)
library(RColorBrewer)



options(stringsAsFactors = FALSE)

source(here("scripts","00_pretty_plots.R"))

## color scheme
cols_manual<-c(  "pink",
                 "#3D0AF2","#6521C1","#67D364","#367C34",
                 "#B921C1","#308AC8",
                 "#C86E30","#C12134","#238443","#810f7c" ,"#02818a","#8c510a" ,"#78c679","#67a9cf",
                 "#3690c0","#8073ac","#b2abd2","#238443","#dd3497","#1d91c0","#addd8e","#014636",
                 "#253494","#081d58","#bf812d","#7a0177" ,"#ce1256","#006d2c" ,
                 "#a50f15","#542788" ,"#e31a1c","#e08214","#ef6548","#fdd49e" ,"#f46d43","#4575b4" ,
                 "#016c59", "#bf812d","#c51b7d","#fb9a99","#fb6a4a","lightgrey","black","#c6dbef","#c7e9c0","pink")


names(cols_manual) <- c("Unidentified_only_UC",
                        "BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine",
                        "Goblet cell", "Paneth cell","Tuft",
                        "Activated B cell","activated DC" ,"Activated T","Arterial endothelial cell" ,"B cell","CD4 T cell",
                        "CD8 T cell","cDC1","cDC2","Cycling B cell","Cycling myeloid cells","Cycling plasma cell","FCER2 B cell","gd T/NK cell",
                        "IgA plasma cell","IgG plasma cell","Lymphatic endothelial cell","Macrophage" ,"mast cells","Memory B cell" ,
                        "Monocyte","pDC" ,"pericyte","S1 fibroblasts","S2 fibroblasts","S4 fibroblasts" ,"TA","Tfh" ,
                        "Treg", "Venous endothelial cell","Glial cell","myofibroblast","Mast cell","unknown_filtered","Doublet","Neonatal CD4 T cell","Neonatal B cell","Neonatal Epithelial")
fillscale_cols_manual <- scale_fill_manual(name="Cell Type",values = cols_manual, drop = T)
colcale_cols_manual <- scale_color_manual(name="Cell Type",values = cols_manual, drop = T)





## this data is filtered genes with expression in less than 3 cells, cells <200 or > 6000 n_feature, percent MT >20 and doublets
# but not normalized or scaled
d10x.primary<-readRDS(here("data","d10x_primary_raw_merged.rds"))

## add cell type labels from split analysis
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

## remove doublets
d10x.primary<-subset(d10x.primary, subset = cluster_ID != "Doublet")
d10x.primary.cd.ctrl<-subset(d10x.primary, subset = orig.ident %in% c("CD","Ctrl"))

table(d10x.primary.cd.ctrl$cluster_ID)
            # 
            # d10x.primary.cd.ctrl_subset<- subset(d10x.primary.cd.ctrl, subset = cluster_ID %in% c("enterocyte","Activated B cell"))
            # 
            # 
            # #################
            # ## Combined analysis with 1 0 in meta data
            # ################
            # write.table(as.matrix(GetAssayData(object = d10x.primary.cd.ctrl_subset, slot = "counts")), 
            #             '../../data/cellphonedb/primary_CDctrl_counts.txt', 
            #             sep = '\t', row.names = T, col.names = T, quote = F)
            # 
            # meta_data<-d10x.primary.cd.ctrl_subset@meta.data[,c("orig.ident","cluster_ID")]
            # meta_data$orig.ident<-as.character(meta_data$orig.ident)
            # meta_data$cell_diagnosis<-paste(meta_data$cluster_ID, meta_data$orig.ident)
            # 
            # meta_data$Cell<-rownames(meta_data)
            # 
            # write.table(meta_data[,c("Cell","cell_diagnosis")], 
            #             '../../data/cellphonedb/primary_CDctrl_meta.txt', 
            #             sep = '\t', row.names = T, col.names = T, quote = F)
            # 
            # 
            # #cellphonedb method statistical_analysis primary_CDctrl_meta.txt primary_CDctrl_counts.txt --counts-data=gene_name
            # #cellphonedb plot dot_plot
            # 

          # 
          # ############################
          # ## CD ctrl separate
          # ############################
          # 
          # #We recommend using normalized count data. Importantly, the user needs to specify whether the data was log-transformed when using the subsampling option
          # d10x.primary.cd.ctrl_subset <- NormalizeData(d10x.primary.cd.ctrl_subset,scale.factor = 10000, normalization.method = "RC")
          # 
          # 
          # d10x.primary.cd_subset<- subset(d10x.primary.cd.ctrl_subset, subset = orig.ident == "CD")
          # d10x.primary.ctrl_subset<- subset(d10x.primary.cd.ctrl_subset, subset = orig.ident == "Ctrl")
          # 
          # 
          # write.table(as.matrix(GetAssayData(object = d10x.primary.cd_subset, slot = "counts")), 
          #             '../../data/cellphonedb/CD_ctrl_seperate/CD/primary_CD_counts.txt', 
          #             sep = '\t', row.names = T, col.names = T, quote = F)
          # write.table(as.matrix(GetAssayData(object = d10x.primary.ctrl_subset, slot = "counts")), 
          #             '../../data/cellphonedb/CD_ctrl_seperate/CTRL/primary_CTRL_counts.txt', 
          #             sep = '\t', row.names = T, col.names = T, quote = F)
          # 
          # meta_data_cd<-d10x.primary.cd_subset@meta.data[,c("orig.ident","cluster_ID")]
          # meta_data_cd$orig.ident<-as.character(meta_data_cd$orig.ident)
          # meta_data_cd$Cell<-rownames(meta_data_cd)
          # 
          # write.table(meta_data_cd[,c("Cell","cluster_ID")], 
          #             '../../data/cellphonedb/CD_ctrl_seperate/CD/primary_CD_meta.txt', 
          #             sep = '\t', row.names = T, col.names = T, quote = F)
          # 
          # meta_data_ctrl<-d10x.primary.ctrl_subset@meta.data[,c("orig.ident","cluster_ID")]
          # meta_data_ctrl$orig.ident<-as.character(meta_data_ctrl$orig.ident)
          # meta_data_ctrl$Cell<-rownames(meta_data_ctrl)
          # 
          # write.table(meta_data_ctrl[,c("Cell","cluster_ID")], 
          #             '../../data/cellphonedb/CD_ctrl_seperate/CTRL/primary_CTRL_meta.txt', 
          #             sep = '\t', row.names = T, col.names = T, quote = F)
          # 
          
          #cellphonedb method statistical_analysis primary_CD_meta.txt primary_CD_counts.txt --counts-data=gene_name
          #cellphonedb plot dot_plot
          
          #cellphonedb method statistical_analysis primary_CTRL_meta.txt primary_CTRL_counts.txt --counts-data=gene_name
          #cellphonedb plot dot_plot


            # 
            # 
            # CDmeans<-read.table(file=here("../../data/cellphonedb/CD_ctrl_seperate/CD/out/means.txt"), sep="\t", header=T)
            # CDpvals<-read.table(file=here("../../data/cellphonedb/CD_ctrl_seperate/CD/out/pvalues.txt"), sep="\t", header=T)
            # CDsigmeans<-read.table(file=here("../../data/cellphonedb/CD_ctrl_seperate/CD/out/significant_means.txt"), sep="\t", header=T)
            # 
            # CTRLsigmeans<-read.table(file=here("../../data/cellphonedb/CD_ctrl_seperate/CTRL/out/significant_means.txt"), sep="\t", header=T)
            # CTRLpvals<-read.table(file=here("../../data/cellphonedb/CD_ctrl_seperate/CTRL/out/pvalues.txt"), sep="\t", header=T)
            # 
            # 
            # ## confrim how means done
            # cdcnt<-as.matrix(GetAssayData(object = d10x.primary.cd_subset, slot = "counts"))
            # cdcnt[grep("CADM1|NECTIN3", rownames(cdcnt)),]
            # 
            # tapply(cdcnt[which(rownames(cdcnt)=="CADM1"),], meta_data_cd$cluster_ID, mean)
            # tapply(cdcnt[which(rownames(cdcnt)=="NECTIN3"),], meta_data_cd$cluster_ID, mean)
            # 
            # 
            # ## CTRL CD differences
            # nrow(CDsigmeans)
            # CDsigmeans_sig<-CDsigmeans[which(sapply(1:nrow(CDsigmeans), function(x) sum(!(is.na(CDsigmeans[x,13:16]))))>0),]
            # nrow(CDsigmeans_sig)
            # 
            # nrow(CTRLsigmeans)
            # CTRLsigmeans_sig<-CTRLsigmeans[which(sapply(1:nrow(CTRLsigmeans), function(x) sum(!(is.na(CTRLsigmeans[x,13:16]))))>0),]
            # nrow(CTRLsigmeans_sig)
            # 
            # length(intersect(CTRLsigmeans_sig$interacting_pair, CDsigmeans_sig$interacting_pair))
            # 
            # CDsigmeans_unique<-CDsigmeans_sig[which(!(CDsigmeans_sig$interacting_pair%in%CTRLsigmeans_sig$interacting_pair)),]
            # 
            # CDpvals_unique<-CDpvals[which(!(CDpvals$interacting_pair%in%CTRLpvals$interacting_pair)),]
            # 
            # 
            # ## why are some pairs not present in both files
            # # Are they many not detected in one group?
            # # Yes seems like they are excluded from P value if 0 in all cells
            # not_in_ctrl<-CDsigmeans[which(!(CDsigmeans$interacting_pair%in%CTRLsigmeans$interacting_pair)),]
            # #CCL3     CCR1
            # ctrlcnt<-as.matrix(GetAssayData(object = d10x.primary.ctrl_subset, slot = "counts"))
            # ctrlcnt[grep("CCL3|CCR1", rownames(ctrlcnt)),1:5]
            # ctrlcnt[grep("CADM1|NECTIN3", rownames(ctrlcnt)),1:5]
            # 
            # tapply(ctrlcnt[which(rownames(cdcnt)=="CCL3"),], meta_data_ctrl$cluster_ID, mean)
            # tapply(ctrlcnt[which(rownames(cdcnt)=="CCR1"),], meta_data_ctrl$cluster_ID, mean)
            # 
            # tapply(ctrlcnt[which(rownames(cdcnt)=="TNFSF14"),], meta_data_ctrl$cluster_ID, mean)
            # tapply(ctrlcnt[which(rownames(cdcnt)=="TNFRSF14"),], meta_data_ctrl$cluster_ID, mean)
            # 
            # tapply(ctrlcnt[which(rownames(cdcnt)=="CXCL11"),], meta_data_ctrl$cluster_ID, mean)
            # tapply(ctrlcnt[which(rownames(cdcnt)=="DPP4"),], meta_data_ctrl$cluster_ID, mean)
            # 
            # 



#############################################
## Full analysis
#############################################
d10x.primary.cd.ctrl_subset_epi_Tcell<- subset(d10x.primary.cd.ctrl, subset = cluster_ID %in% c("BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine", "Goblet cell", "Paneth cell","Tuft","TA",
                                                                                      "Activated T","CD4 T cell","CD8 T cell","gd T/NK cell","Tfh","Treg"))

## CD ctrl separate

#We recommend using normalized count data. Importantly, the user needs to specify whether the data was log-transformed when using the subsampling option
d10x.primary.cd.ctrl_subset_epi_Tcell <- NormalizeData(d10x.primary.cd.ctrl_subset_epi_Tcell,scale.factor = 10000, normalization.method = "RC")


d10x.primary.cd_subset<- subset(d10x.primary.cd.ctrl_subset_epi_Tcell, subset = orig.ident == "CD")
d10x.primary.ctrl_subset<- subset(d10x.primary.cd.ctrl_subset_epi_Tcell, subset = orig.ident == "Ctrl")


write.table(as.matrix(GetAssayData(object = d10x.primary.cd_subset, slot = "counts")), 
            file=paste(here("data/"),"primary_CD_counts.txt", sep=""), 
            sep = '\t', row.names = T, col.names = T, quote = F)
write.table(as.matrix(GetAssayData(object = d10x.primary.ctrl_subset, slot = "counts")), 
            file=paste(here("data/"),"primary_CTRL_counts.txt", sep=""), 
            sep = '\t', row.names = T, col.names = T, quote = F)

meta_data_cd<-d10x.primary.cd_subset@meta.data[,c("orig.ident","cluster_ID")]
meta_data_cd$orig.ident<-as.character(meta_data_cd$orig.ident)
meta_data_cd$Cell<-rownames(meta_data_cd)

write.table(meta_data_cd[,c("Cell","cluster_ID")], 
            file=paste(here("data/"),"primary_CD_meta.txt", sep=""),  
            sep = '\t', row.names = T, col.names = T, quote = F)

meta_data_ctrl<-d10x.primary.ctrl_subset@meta.data[,c("orig.ident","cluster_ID")]
meta_data_ctrl$orig.ident<-as.character(meta_data_ctrl$orig.ident)
meta_data_ctrl$Cell<-rownames(meta_data_ctrl)

write.table(meta_data_ctrl[,c("Cell","cluster_ID")], 
            file=paste(here("data/"),"primary_CTRL_meta.txt", sep=""), 
            sep = '\t', row.names = T, col.names = T, quote = F)

#scp /home/redgar/Documents/ibd/data/cellphonedb/CD_ctrl_seperate/CTRL/* redgar@ebi:/nfs/research1/zerbino/redgar/ibd/data/cellphonedb/CD_ctrl_seperate/CTRL
#scp /home/redgar/Documents/ibd/data/cellphonedb/CD_ctrl_seperate/CD/* redgar@ebi:/nfs/research1/zerbino/redgar/ibd/data/cellphonedb/CD_ctrl_seperate/CD

# also in cellphonedb_primary_CD.sh
#cellphonedb method statistical_analysis primary_CD_meta.txt primary_CD_counts.txt --counts-data=gene_name 
#cellphonedb plot dot_plot

## also in cellphonedb_primary_CTRL.sh
#cellphonedb method statistical_analysis primary_CTRL_meta.txt primary_CTRL_counts.txt --counts-data=gene_name
#cellphonedb plot dot_plot

  
################################################################################################
###########
## Full Results
###########
CDmeans<-read.table(file=here("data/cellphone_db/CD/out/means.txt"), sep="\t", header=T)
CDpvals<-read.table(file=here("data/cellphone_db/CD/out/pvalues.txt"), sep="\t", header=T)
CDsigmeans<-read.table(file=here("data/cellphone_db/CD/out/significant_means.txt"), sep="\t", header=T)

CTRLmeans<-read.table(file=here("data/cellphone_db/CTRL/out/means.txt"), sep="\t", header=T)
CTRLsigmeans<-read.table(file=here("data/cellphone_db/CTRL/out/significant_means.txt"), sep="\t", header=T)
CTRLpvals<-read.table(file=here("data/cellphone_db/CTRL/out/pvalues.txt"), sep="\t", header=T)



## CTRL CD differences
nrow(CDsigmeans)
CDsigmeans_sig<-CDsigmeans[which(sapply(1:nrow(CDsigmeans), function(x) sum(!(is.na(CDsigmeans[x,13:ncol(CDsigmeans)]))))>0),]
nrow(CDsigmeans_sig)
CDsigmeans_sig$interacting_pair

nrow(CTRLsigmeans)
CTRLsigmeans_sig<-CTRLsigmeans[which(sapply(1:nrow(CTRLsigmeans), function(x) sum(!(is.na(CTRLsigmeans[x,13:ncol(CDsigmeans)]))))>0),]
nrow(CTRLsigmeans_sig)
CTRLsigmeans_sig$interacting_pair

length(intersect(CTRLsigmeans_sig$interacting_pair, CDsigmeans_sig$interacting_pair))

CDsigmeans_unique<-CDsigmeans_sig[which(!(CDsigmeans_sig$interacting_pair%in%CTRLsigmeans_sig$interacting_pair)),]
CDsigmeans_unique$interacting_pair
CDsigmeans_unique[,1:15]

lapply(1:nrow(CDsigmeans_unique), function(x){
  paste(CDsigmeans_unique$interacting_pair[x],":",paste0(colnames(CDsigmeans_unique)[which(!(is.na(CDsigmeans_unique[x,])))][-c(1:12)], collapse = ", "))
})



## fixing colnames

cell_types_spaces<-c("BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine", "Goblet cell", "Paneth cell","Tuft","TA","Activated T","CD4 T cell","CD8 T cell","gd T/NK cell","Tfh","Treg")


combos_spaces<-combinations(n = 15, r = 2, v = cell_types_spaces, repeats.allowed = T)
easytouse_names1<-apply(combos_spaces,1, function(pair) paste(pair, collapse="."))
easytouse_names2<-apply(combos_spaces[,c(2,1)],1, function(pair) paste(pair, collapse="."))
easytouse_names<-c(easytouse_names1,easytouse_names2)

cell_types<-gsub(" ",".",cell_types_spaces)
combos<-combinations(n = 15, r = 2, v = cell_types, repeats.allowed = T)
hardtouse_names1<-apply(combos,1, function(pair) paste(pair, collapse="."))
hardtouse_names2<-apply(combos[,c(2,1)],1, function(pair) paste(pair, collapse="."))
hardtouse_names<-c(hardtouse_names1,hardtouse_names2)
hardtouse_names<-gsub("/",".",hardtouse_names)

length(colnames(CDsigmeans_unique)[-c(1:12)])
length(which(colnames(CDsigmeans_unique)[-c(1:12)]%in%hardtouse_names))
colnames(CDsigmeans_unique)[-c(1:12)][which(!(colnames(CDsigmeans_unique)[-c(1:12)]%in%hardtouse_names))]


#change colnames
change<-data.frame(easy=easytouse_names, hard=hardtouse_names)
change<-change[match(colnames(CDsigmeans_unique)[-c(1:12)],change$hard),]

colnames(CDsigmeans_unique)[-c(1:12)]<-change$easy

# interactions are not symmetric. Partner A expression is considered on the first cluster, and partner B expression is considered on the second cluster. In other words:
melted_interactions<-do.call(rbind, lapply(1:nrow(CDsigmeans_unique), function(x){
  melted<-as.data.frame(do.call(rbind,strsplit(colnames(CDsigmeans_unique)[which(!(is.na(CDsigmeans_unique[x,])))][-c(1:12)], "[.]")))
  colnames(melted)<-c("cell_a","cell_b")
  melted$gene_a<-CDsigmeans_unique$gene_a[x]
  melted$gene_b<-CDsigmeans_unique$gene_b[x]
  melted[,c(3,4,1,2)]
}))


################
## Heat mapepi immune
################

## epi immune
CD_interactions_sigdiff<-CDsigmeans_unique$interacting_pair


index_cell_cell<-cbind(change, as.data.frame(do.call(rbind,strsplit(change$easy, "[.]"))))
colnames(index_cell_cell)<-c("easy", "hard", "cell_a","cell_b")


immune<-c("Activated T","CD4 T cell","CD8 T cell","gd T/NK cell","Tfh","Treg")
epi<-c("BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine", "Goblet cell", "Paneth cell","Tuft","TA")



epi_immune<-do.call(rbind, lapply(1:nrow(index_cell_cell), function(x){
  if((index_cell_cell$cell_a[x]%in%immune & index_cell_cell$cell_b[x]%in%epi) | (index_cell_cell$cell_a[x]%in%epi & index_cell_cell$cell_b[x]%in%immune)){
    index_cell_cell[x,]
  }
}))

cell_cell_to_pull<-epi_immune$hard

plt_CD_mean<-CDmeans[which(CDmeans$interacting_pair%in%CD_interactions_sigdiff), which(colnames(CDmeans)%in%c("interacting_pair",cell_cell_to_pull))]
plt_CD_pval<-CDpvals[which(CDpvals$interacting_pair%in%CD_interactions_sigdiff), which(colnames(CDpvals)%in%c("interacting_pair",cell_cell_to_pull))]

plt_CD_mean<-melt(plt_CD_mean)
colnames(plt_CD_mean)<-c("interacting_pair","cells","mean")
plt_CD_pval<-melt(plt_CD_pval)
colnames(plt_CD_pval)<-c("interacting_pair","cells","pval")

plt_CD<-merge(plt_CD_mean,plt_CD_pval, by=c("cells","interacting_pair"))

## sig between epi and immune
epi_immune_pairs<-unique(plt_CD$interacting_pair[which(plt_CD$pval<0.05)])
plt_CD<-plt_CD[which(plt_CD$interacting_pair%in%epi_immune_pairs),]

## cell type label
index_cell_cell_plt<-index_cell_cell[which(index_cell_cell$hard%in%plt_CD$cells),]
index_cell_cell_plt<-index_cell_cell_plt[match(plt_CD$cells,index_cell_cell_plt$hard),]

a<-index_cell_cell_plt[,c("hard","cell_a")]
colnames(a)<-c("hard","cell")
a$posistion<-"1"
b<-index_cell_cell_plt[,c("hard","cell_b")]
colnames(b)<-c("hard","cell")
b$posistion<-"2"

index_cell_cell_plt_tidy<-rbind(a,b)

## cell order
index_cell_cell_plt$cell_a<-factor(index_cell_cell_plt$cell_a, levels = cell_types_spaces)
index_cell_cell_plt$cell_b<-factor(index_cell_cell_plt$cell_b, levels = cell_types_spaces)

index_cell_cell_plt_order<-index_cell_cell_plt[order(index_cell_cell_plt$cell_a, index_cell_cell_plt$cell_b),]
index_cell_cell_plt_order_hard<-unique(index_cell_cell_plt_order$hard)

## cell label legend
index_cell_cell_plt_tidy$hard<-factor(index_cell_cell_plt_tidy$hard, levels=index_cell_cell_plt_order_hard)
rectanle<-ggplot(index_cell_cell_plt_tidy, aes(hard,posistion))+geom_tile(aes(fill=cell))+theme(legend.position = "bottom",
                                                                                 axis.title=element_blank(),
                                                                                 axis.text=element_blank(),
                                                                                 axis.ticks=element_blank(),
                                                                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                 panel.background = element_blank(), axis.line = element_blank(),
                                                                                 panel.border = element_blank())+fillscale_cols_manual



## mean dot plot
plt_CD$cells<-factor(plt_CD$cells, levels=index_cell_cell_plt_order_hard)

plt_CD$pval_binary<-"Not Significant"
plt_CD$pval_binary[which(plt_CD$pval<0.05)]<-"Significant"

#myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

myPalette <- colorRampPalette(colorspace::sequential_hcl(11, h = c(300, 75), c = c(35, 95), l = c(15, 90), power = c(0.8, 1.2)))

plt_CD$diagnosis<-"CD"

              dotplt<-ggplot(plt_CD, aes(cells,interacting_pair))+geom_point(aes(fill=log(mean), size=pval_binary), shape=21, color="grey90")+theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
                scale_size_manual(values=c(2,6), name="Interaction\nsignificance")+geom_vline(xintercept=54.5, color="grey")+scale_fill_gradientn(colours = myPalette(100),  name="Interacting pair\nmean expression\n(log)", na.value = "white")


              plot_grid(dotplt,rectanle,
                        ncol = 1, align="v", axis = "lr",
                        rel_heights = c(4,1))







##### controls

plt_CTRL_mean<-CTRLmeans[which(CTRLmeans$interacting_pair%in%CD_interactions_sigdiff), which(colnames(CTRLmeans)%in%c("interacting_pair",cell_cell_to_pull))]
plt_CTRL_pval<-CTRLpvals[which(CTRLpvals$interacting_pair%in%CD_interactions_sigdiff), which(colnames(CTRLpvals)%in%c("interacting_pair",cell_cell_to_pull))]

plt_CTRL_mean<-melt(plt_CTRL_mean)
colnames(plt_CTRL_mean)<-c("interacting_pair","cells","mean")
plt_CTRL_pval<-melt(plt_CTRL_pval)
colnames(plt_CTRL_pval)<-c("interacting_pair","cells","pval")

plt_CTRL<-merge(plt_CTRL_mean,plt_CTRL_pval, by=c("cells","interacting_pair"))

## sig between epi and immune
plt_CTRL<-plt_CTRL[which(plt_CTRL$interacting_pair%in%epi_immune_pairs),]




## mean dot plot
plt_CTRL$cells<-factor(plt_CTRL$cells, levels=index_cell_cell_plt_order_hard)

plt_CTRL$pval_binary<-"Not Significant"
plt_CTRL$pval_binary[which(plt_CTRL$pval<0.05)]<-"Significant"

plt_CTRL$diagnosis<-"Control"

plt_diagnosis<-rbind(plt_CD, plt_CTRL)
              #
              # dotplt<-ggplot(plt_diagnosis, aes(cells,interacting_pair))+geom_point(aes(fill=log(mean), size=pval_binary), shape=21)+theme_bw()+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
              #   scale_size_manual(values=c(2,6), name="Interaction\nsignificance")+geom_vline(xintercept=54.5, color="grey")+scale_fill_gradientn(colours = myPalette(100),  name="Interacting pair\nmean expression\n(log)", na.value = "white")


              # dotplt<-ggplot(plt_diagnosis, aes(cells,diagnosis))+geom_point(aes(fill=log(mean), size=pval_binary), shape=21)+facet_grid(interacting_pair~.)+
              #   theme_bw()+   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text.y = element_text(angle = 0))+
              #   scale_size_manual(values=c(2,6), name="Interaction\nsignificance")+geom_vline(xintercept=54.5, color="grey")+scale_fill_gradientn(colours = myPalette(100),  name="Interacting pair\nmean expression\n(log)", na.value = "white")
              #
              #
              # plot_grid(dotplt,rectanle,
              #           ncol = 1, align="v", axis = "lr",
              #           rel_heights = c(4,1))
              #



## split things up a bit


plot_some_interactions<-function(inter_of_interest){
  plt_diagnosis<-rbind(plt_CD, plt_CTRL)
  plt_diagnosis<-plt_diagnosis[which(plt_diagnosis$interacting_pair%in%inter_of_interest),]

  dotplt<-ggplot(plt_diagnosis, aes(cells,diagnosis))+geom_point(aes(fill=log(mean), size=pval_binary), shape=21, color="grey90")+facet_grid(interacting_pair~.,switch="y")+th+
    theme_bw()+   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text.y = element_text(angle = 180),panel.spacing = unit(0, "lines"))+scale_y_discrete(position = "right")+xlab("Cell Pair")+
    scale_size_manual(values=c(2,6), name="Interaction\nsignificance")+geom_vline(xintercept=54.5, color="grey", size=0.5)+scale_fill_gradientn(colours = myPalette(100),  name="Interacting pair\nmean expression\n(log)", na.value = "white")


  plot_grid(dotplt,rectanle,
            ncol = 1, align="v", axis = "lr",
            rel_heights = c(4,1))}

#all
plot_some_interactions(unique(plt_CD$interacting_pair))
ggsave(file=here("../../../codon/scRNAseq_codon/figs/jpeg/all_CD_CTRL_epiimmune.jpeg"), width=25, height=25)
ggsave(file=here("../../../codon/scRNAseq_codon/figs/all_CD_CTRL_epiimmune.pdf"), width=25, height=25)

colagen<-unique(plt_CD$interacting_pair[grep("COL",plt_CD$interacting_pair)])
plot_some_interactions(colagen)
ggsave(file=here("../../../codon/scRNAseq_codon/figs/jpeg/all_CD_CTRL_epiimmune_COL.jpeg"), width=25, height=10)
ggsave(file=here("../../../codon/scRNAseq_codon/figs/all_CD_CTRL_epiimmune_COL.pdf"), width=25, height=10)


ICAM<-unique(plt_CD$interacting_pair[grep("ICAM",plt_CD$interacting_pair)])
plot_some_interactions(ICAM)
ggsave(file=here("../../../codon/scRNAseq_codon/figs/jpeg/all_CD_CTRL_epiimmune_ICAM.jpeg"), width=25, height=8)
ggsave(file=here("../../../codon/scRNAseq_codon/figs/all_CD_CTRL_epiimmune_ICAM.pdf"), width=25, height=8)


CXCL<-unique(plt_CD$interacting_pair[grep("CXCL",plt_CD$interacting_pair)])
plot_some_interactions(CXCL)
ggsave(file=here("../../../codon/scRNAseq_codon/figs/jpeg/all_CD_CTRL_epiimmune_CXCL.jpeg"), width=25, height=6)
ggsave(file=here("../../../codon/scRNAseq_codon/figs/all_CD_CTRL_epiimmune_CXCL.pdf"), width=25, height=6)



other<-unique(plt_CD$interacting_pair[-grep("CXCL|ICAM|COL",plt_CD$interacting_pair)])
plot_some_interactions(other)
ggsave(file=here("../../../codon/scRNAseq_codon/figs/jpeg/all_CD_CTRL_epiimmune_other.jpeg"), width=25, height=15)
ggsave(file=here("../../../codon/scRNAseq_codon/figs/all_CD_CTRL_epiimmune_other.pdf"), width=25, height=15)



################
## Heat map epi epi
################
# 
# CD_interactions_sigdiff<-CDsigmeans_unique$interacting_pair
# 
# index_cell_cell<-cbind(change, as.data.frame(do.call(rbind,strsplit(change$easy, "[.]"))))
# colnames(index_cell_cell)<-c("easy", "hard", "cell_a","cell_b")
# 
# epi<-c("BEST4 enterocyte","crypt","early enterocyte","enterocyte","enteroendocrine", "Goblet cell", "Paneth cell","Tuft","TA")
# 
# epi_epi<-do.call(rbind, lapply(1:nrow(index_cell_cell), function(x){
#   if((index_cell_cell$cell_a[x]%in%epi & index_cell_cell$cell_b[x]%in%epi) | (index_cell_cell$cell_a[x]%in%epi & index_cell_cell$cell_b[x]%in%epi)){
#     index_cell_cell[x,]
#   }
# }))
# 
# cell_cell_to_pull<-epi_epi$hard
# 
# plt_CD_mean<-CDmeans[which(CDmeans$interacting_pair%in%CD_interactions_sigdiff), which(colnames(CDmeans)%in%c("interacting_pair",cell_cell_to_pull))]
# plt_CD_pval<-CDpvals[which(CDpvals$interacting_pair%in%CD_interactions_sigdiff), which(colnames(CDpvals)%in%c("interacting_pair",cell_cell_to_pull))]
# 
# plt_CD_mean<-melt(plt_CD_mean)
# colnames(plt_CD_mean)<-c("interacting_pair","cells","mean")
# plt_CD_pval<-melt(plt_CD_pval)
# colnames(plt_CD_pval)<-c("interacting_pair","cells","pval")
# 
# plt_CD<-merge(plt_CD_mean,plt_CD_pval, by=c("cells","interacting_pair"))
# 
# ## sig between epi and immune
# epi_epi_pairs<-unique(plt_CD$interacting_pair[which(plt_CD$pval<0.05)])
# plt_CD<-plt_CD[which(plt_CD$interacting_pair%in%epi_epi_pairs),]
# 
# ## cell type label
# index_cell_cell_plt<-index_cell_cell[which(index_cell_cell$hard%in%plt_CD$cells),]
# index_cell_cell_plt<-index_cell_cell_plt[match(plt_CD$cells,index_cell_cell_plt$hard),]
# 
# a<-index_cell_cell_plt[,c("hard","cell_a")]
# colnames(a)<-c("hard","cell")
# a$posistion<-"1"
# b<-index_cell_cell_plt[,c("hard","cell_b")]
# colnames(b)<-c("hard","cell")
# b$posistion<-"2"
# 
# index_cell_cell_plt_tidy<-rbind(a,b)
# 
# ## cell order
# index_cell_cell_plt$cell_a<-factor(index_cell_cell_plt$cell_a, levels = cell_types_spaces)
# index_cell_cell_plt$cell_b<-factor(index_cell_cell_plt$cell_b, levels = cell_types_spaces)
# 
# index_cell_cell_plt_order<-index_cell_cell_plt[order(index_cell_cell_plt$cell_a, index_cell_cell_plt$cell_b),]
# index_cell_cell_plt_order_hard<-unique(index_cell_cell_plt_order$hard)
# 
# ## cell label legend
# index_cell_cell_plt_tidy$hard<-factor(index_cell_cell_plt_tidy$hard, levels=index_cell_cell_plt_order_hard)
# rectanle<-ggplot(index_cell_cell_plt_tidy, aes(hard,posistion))+geom_tile(aes(fill=cell))+theme(legend.position = "bottom",
#                                                                                                 axis.title=element_blank(),
#                                                                                                 axis.text=element_blank(),
#                                                                                                 axis.ticks=element_blank(),
#                                                                                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                                                                                                 panel.background = element_blank(), axis.line = element_blank(),
#                                                                                                 panel.border = element_blank())+fillscale_cols_manual
# 
# 
# 
# ## mean dot plot
# plt_CD$cells<-factor(plt_CD$cells, levels=index_cell_cell_plt_order_hard)
# 
# plt_CD$pval_binary<-"Not Significant"
# plt_CD$pval_binary[which(plt_CD$pval<0.05)]<-"Significant"
# 
# #myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
# 
# myPalette <- colorRampPalette(colorspace::sequential_hcl(11, h = c(300, 75), c = c(35, 95), l = c(15, 90), power = c(0.8, 1.2)))
# 
# plt_CD$diagnosis<-"CD"
# 
# 
# ##### controls
# plt_CTRL_mean<-CTRLmeans[which(CTRLmeans$interacting_pair%in%CD_interactions_sigdiff), which(colnames(CTRLmeans)%in%c("interacting_pair",cell_cell_to_pull))]
# plt_CTRL_pval<-CTRLpvals[which(CTRLpvals$interacting_pair%in%CD_interactions_sigdiff), which(colnames(CTRLpvals)%in%c("interacting_pair",cell_cell_to_pull))]
# 
# plt_CTRL_mean<-melt(plt_CTRL_mean)
# colnames(plt_CTRL_mean)<-c("interacting_pair","cells","mean")
# plt_CTRL_pval<-melt(plt_CTRL_pval)
# colnames(plt_CTRL_pval)<-c("interacting_pair","cells","pval")
# 
# plt_CTRL<-merge(plt_CTRL_mean,plt_CTRL_pval, by=c("cells","interacting_pair"))
# 
# ## sig between epi and immune
# plt_CTRL<-plt_CTRL[which(plt_CTRL$interacting_pair%in%epi_epi_pairs),]
# 
# 
# 
# 
# ## mean dot plot
# plt_CTRL$cells<-factor(plt_CTRL$cells, levels=index_cell_cell_plt_order_hard)
# 
# plt_CTRL$pval_binary<-"Not Significant"
# plt_CTRL$pval_binary[which(plt_CTRL$pval<0.05)]<-"Significant"
# 
# plt_CTRL$diagnosis<-"Control"
# 
# plt_diagnosis<-rbind(plt_CD, plt_CTRL)
# 
# 
# ## split things up a bit
# plot_some_interactions<-function(inter_of_interest){
#   plt_diagnosis<-rbind(plt_CD, plt_CTRL)
#   plt_diagnosis<-plt_diagnosis[which(plt_diagnosis$interacting_pair%in%inter_of_interest),]
# 
#   dotplt<-ggplot(plt_diagnosis, aes(cells,diagnosis))+geom_point(aes(fill=log(mean), size=pval_binary), shape=21, color="grey90")+facet_grid(interacting_pair~.,switch="y")+
#     theme_bw()+   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),strip.text.y = element_text(angle = 180),panel.spacing = unit(0, "lines"))+scale_y_discrete(position = "right")+xlab("Cell Pair")+
#     scale_size_manual(values=c(2,6), name="Interaction\nsignificance")+scale_fill_gradientn(colours = myPalette(100),  name="Interacting pair\nmean expression\n(log)", na.value = "white")
# 
# 
#   plot_grid(dotplt,rectanle,
#             ncol = 1, align="v", axis = "lr",
#             rel_heights = c(4,1))}
# 
# #all
# plot_some_interactions(unique(plt_CD$interacting_pair))
# 
# colagen<-unique(plt_CD$interacting_pair[grep("COL",plt_CD$interacting_pair)])
# plot_some_interactions(colagen)
# 
# FGF<-unique(plt_CD$interacting_pair[grep("FGF",plt_CD$interacting_pair)])
# plot_some_interactions(FGF)
# 
# CXCL<-unique(plt_CD$interacting_pair[grep("CXCL",plt_CD$interacting_pair)])
# plot_some_interactions(CXCL)
# 
# nothcwnt<-unique(plt_CD$interacting_pair[grep("NOTCH|WNT",plt_CD$interacting_pair)])
# plot_some_interactions(nothcwnt)
# 
# ack_ccl<-unique(plt_CD$interacting_pair[grep("ACK|CCL",plt_CD$interacting_pair)])
# plot_some_interactions(ack_ccl)
# 
# bmp<-unique(plt_CD$interacting_pair[grep("BMP",plt_CD$interacting_pair)])
# plot_some_interactions(bmp)
# 
# other<-unique(plt_CD$interacting_pair[-grep("CXCL|FGF|COL|NOTCH|WNT|ACK|CCL|BMP",plt_CD$interacting_pair)])
# plot_some_interactions(other)



##################
## does CD have more interactions between epi and immune?
##################
length(CTRLsigmeans_sig$interacting_pair)
length(CDsigmeans_sig$interacting_pair)


CDsigmeans_sig_melt<-melt(CDsigmeans_sig[, which(colnames(CDsigmeans_sig)%in%c("interacting_pair",index_cell_cell$hard))], id="interacting_pair")
colnames(CDsigmeans_sig_melt)<-c("interacting_pair","cells","mean")
CDsigmeans_sig_melt<-CDsigmeans_sig_melt[!is.na(CDsigmeans_sig_melt$mean),]

index_cell_cell$hard2<-gsub(" ",".",paste(index_cell_cell$cell_b, index_cell_cell$cell_a, sep=" "))

sum_CDsigmeans<-as.data.frame(table(CDsigmeans_sig_melt$cells))
sum_CDsigmeans<-merge(sum_CDsigmeans,index_cell_cell, by.x="Var1", by.y="hard")

cell_cell_pairs<-unique(index_cell_cell$hard)
sum_CDsigmeans<-do.call(rbind,lapply(1:length(cell_cell_pairs),function(x){
  df<-sum_CDsigmeans[which(sum_CDsigmeans$Var1==cell_cell_pairs[x] | sum_CDsigmeans$hard2==cell_cell_pairs[x]),]
  data.frame(cell_a=df$cell_a[1], cell_b=df$cell_b[1], interaction_count=sum(df$Freq))}))
sum_CDsigmeans<-sum_CDsigmeans[!duplicated(sum_CDsigmeans),]



CTRLsigmeans_sig_melt<-melt(CTRLsigmeans_sig[, which(colnames(CTRLsigmeans_sig)%in%c("interacting_pair",index_cell_cell$hard))], id="interacting_pair")
colnames(CTRLsigmeans_sig_melt)<-c("interacting_pair","cells","mean")
CTRLsigmeans_sig_melt<-CTRLsigmeans_sig_melt[!is.na(CTRLsigmeans_sig_melt$mean),]

index_cell_cell$hard2<-gsub(" ",".",paste(index_cell_cell$cell_b, index_cell_cell$cell_a, sep=" "))

CTRLsigmeans_sig_melt[x,]
sum_CTRLsigmeans<-as.data.frame(table(CTRLsigmeans_sig_melt$cells))
sum_CTRLsigmeans<-merge(sum_CTRLsigmeans,index_cell_cell, by.x="Var1", by.y="hard")

cell_cell_pairs<-unique(index_cell_cell$hard)
sum_CTRLsigmeans<-do.call(rbind,lapply(1:length(cell_cell_pairs),function(x){
  df<-sum_CTRLsigmeans[which(sum_CTRLsigmeans$Var1==cell_cell_pairs[x] | sum_CTRLsigmeans$hard2==cell_cell_pairs[x]),]
  data.frame(cell_a=df$cell_a[1], cell_b=df$cell_b[1], interaction_count=sum(df$Freq))}))
sum_CTRLsigmeans<-sum_CTRLsigmeans[!duplicated(sum_CTRLsigmeans),]


merge(sum_CDsigmeans, sum_CTRLsigmeans, by=c("cell_a","cell_b"))
# 
# write.csv(sum_CDsigmeans, file="../../output/cellphonedb_CD_edges.csv", quote = F)
# write.csv(sum_CTRLsigmeans, file="../../output/cellphonedb_CTRL_edges.csv", quote = F)
# 
# node_CD<-as.data.frame(tapply(sum_CDsigmeans$interaction_count, sum_CDsigmeans$cell_a, sum))
# colnames(node_CD)<-"inter_number"
# node_CD$cell<-rownames(node_CD)
# node_CD$compartment<-"Immune"
# node_CD$compartment[which(node_CD$cell%in%epi)]<-"Epithelial"
# write.csv(node_CD, file="../../output/cellphonedb_CD_node.csv", quote = F)
# 
# node_CTRL<-as.data.frame(tapply(sum_CTRLsigmeans$interaction_count, sum_CTRLsigmeans$cell_a, sum))
# colnames(node_CTRL)<-"inter_number"
# node_CTRL$cell<-rownames(node_CTRL)
# node_CTRL$compartment<-"Immune"
# node_CTRL$compartment[which(node_CTRL$cell%in%epi)]<-"Epithelial"
# write.csv(node_CTRL, file="../../output/cellphonedb_CTRL_node.csv", quote = F)
# 


####################
##interactions between immune and epi count
####################

combined_interactions<-merge(sum_CDsigmeans, sum_CTRLsigmeans, by=c("cell_a","cell_b"))

combined_immune_epi<-combined_interactions[which( (combined_interactions$cell_a%in%immune & combined_interactions$cell_b%in%epi) | (combined_interactions$cell_a%in%epi & combined_interactions$cell_b%in%immune) ),]
colnames(combined_immune_epi)<-c("cell_a","cell_b","CD_count","CTRL_count")

cor(combined_immune_epi$CD_count, combined_immune_epi$CTRL_count)

combined_immune_epi$diff<-combined_immune_epi$CTRL_count-combined_immune_epi$CD_count

combined_immune_epi_label<-combined_immune_epi[which(abs(combined_immune_epi$diff)>15),]

ggplot(combined_immune_epi, aes(CTRL_count,CD_count))+geom_point()+geom_abline(slope=1, intercept=0)+xlim(0,55)+ylim(0,55)


noflip<-combined_immune_epi[which(combined_immune_epi$cell_a%in%immune),]
flip<-combined_immune_epi[which(combined_immune_epi$cell_a%in%epi),c( "cell_b","cell_a", "CD_count","CTRL_count", "diff")]
colnames(flip)<-c("cell_a","cell_b","CD_count","CTRL_count","diff")
combined_immune_epi_heat<-rbind(noflip, flip)

## when both sides in the receptor ligand are significant
combined_immune_epi_heat<-as.data.frame(combined_immune_epi_heat %>% group_by(cell_a, cell_b) %>% summarise(CD_count = sum(CD_count),
                                              CTRL_count = sum(CTRL_count),
                                              diff = sum(diff)))

## top differences
combined_immune_epi_heat[which(abs(combined_immune_epi_heat$diff)>15),]


combined_immune_epi_heat$cell_b<-factor(combined_immune_epi_heat$cell_b,
                                        levels=c("crypt","Tuft","Paneth cell","TA","early enterocyte",
                                                 "enteroendocrine","BEST4 enterocyte","Goblet cell","enterocyte"))

levels(combined_immune_epi_heat$cell_b)<-c("Stem","Tuft","Paneth","TA","early enterocyte",
                                           "Enteroendocrine","BEST4 enterocyte","Goblet","Enterocyte")

ggplot(combined_immune_epi_heat, aes(cell_a,cell_b, fill=diff))+geom_tile(color="black")+
  scale_fill_gradient2(low ="#b2182b", mid = "white", high = "#2166ac", midpoint = 0, name="Interaction\nNumber\n(Control-CD) ")+ 
  geom_text(aes(label=paste(CTRL_count,"-",CD_count)), size=3)+
  theme(axis.text = element_text(size =10, color="black"),
        axis.text.x = element_text(),
       axis.title = element_blank(),
       legend.text = element_text(size =10),
       legend.title = element_text(size =10),
       panel.background = element_blank(),
       panel.grid.major = element_blank(),
       panel.grid.minor = element_blank())+th_present+xlab("")+ylab("")
ggsave(file=here("../../../codon/scRNAseq_codon/figs/jpeg/connection_diff_immune_epi.jpeg"), width=10, height=8)
ggsave(file=here("../../../codon/scRNAseq_codon/figs/connection_diff_immune_epi.pdf"), width=10, height=8)

wilcox.test(combined_immune_epi_heat$CD_count, combined_immune_epi_heat$CTRL_count)
