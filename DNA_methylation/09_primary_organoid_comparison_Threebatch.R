#'---
#'title: Organoid and Primary differenital DNAm comparison
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
#library(limma)
library(here)
library(scales)
library(cowplot)
library(dplyr)


options(stringsAsFactors = FALSE)


#' functions
source(here("general_functions","00_pretty_plots.R"))
source(here("general_functions","00_Heat_scree_plot_generic.R"))
source(here("general_functions","00_DNAm_volcano.R"))


#' ## Load DNAm data
load(file=paste(here("DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))
load(here("DNAm/data/ibd_adjusted_combatted.Rdata"))

# load(file=paste(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))
# load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/ibd_adjusted_combatted.Rdata"))

## fix 202 and 203
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605870061_R03C01")]<-"hold"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605880096_R04C01")]<-"205605870061_R03C01"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="hold")]<-"205605880096_R04C01"




#' ## delta beta correlation plot
manhattan_data<-function(CpG, pvalue, fdr, db, diagnosis, tissue, sample.type){
  data.frame(CpG=CpG, db=db, p.value=pvalue, fdr=fdr, diagnosis=diagnosis, tissue=tissue, sample.type=sample.type)
}

sig_hits<-function(manhattan_data_df){
  manhattan_data_df[which(manhattan_data_df$fdr<0.05 & abs(manhattan_data_df$db)>0.05),]
}

load(file=here("DNAm/data/diff_DNAm_grouped_casectrl.RData"))
#load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/diff_DNAm_grouped_casectrl.RData"))

diff_meth_TI_primary<-diff_meth_TI_grouped
diff_meth_SC_primary<-diff_meth_SC_grouped

manhattan_TI_UC_primary<-manhattan_data(diff_meth_TI_primary$CpG, diff_meth_TI_primary$p.value_UC,diff_meth_TI_primary$adjusted_p_UC, diff_meth_TI_primary$db_ctrl_UC, "UC", "TI", "primary")
manhattan_SC_UC_primary<-manhattan_data(diff_meth_SC_primary$CpG, diff_meth_SC_primary$p.value_UC,diff_meth_SC_primary$adjusted_p_UC, diff_meth_SC_primary$db_ctrl_UC, "UC", "SC", "primary")
manhattan_SC_CD_primary<-manhattan_data(diff_meth_SC_primary$CpG, diff_meth_SC_primary$p.value_CD,diff_meth_SC_primary$adjusted_p_CD, diff_meth_SC_primary$db_ctrl_CD, "CD", "SC", "primary")

load(here("DNAm/data/organoid_diff_DNAm_grouped_casectrl_ThreeBatch.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/organoid_diff_DNAm_grouped_casectrl_ThreeBatch.RData"))

diff_meth_TI_organoid<-diff_meth_TI_grouped
diff_meth_SC_organoid<-diff_meth_SC_grouped

manhattan_TI_UC_organoid<-manhattan_data(diff_meth_TI_organoid$CpG, diff_meth_TI_organoid$p.value_UC, diff_meth_TI_organoid$adjusted_p_UC, diff_meth_TI_organoid$db_ctrl_UC, "UC", "TI", "organoid")
manhattan_SC_UC_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_UC, diff_meth_SC_organoid$adjusted_p_UC, diff_meth_SC_organoid$db_ctrl_UC, "UC", "SC", "organoid")
manhattan_SC_CD_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_CD, diff_meth_SC_organoid$adjusted_p_CD, diff_meth_SC_organoid$db_ctrl_CD, "CD", "SC", "organoid")

load(file=here("DNAm/data/diff_DNAm_TI_UCCTRL.RData"))
#load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/diff_DNAm_TI_UCCTRL.RData"))

load(here("DNAm/data/diff_DNAm_TI_UCCTRL_organoid_ThreeBatch.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/diff_DNAm_TI_UCCTRL_organoid_ThreeBatch.RData"))

diff_meth_TI_UCCTRL$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL))
diff_meth_TI_UCCTRL_organoid$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL_organoid$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL_organoid))

manhattan_TI_CD_primary<-manhattan_data(rownames(diff_meth_TI_UCCTRL), diff_meth_TI_UCCTRL$p.value, diff_meth_TI_UCCTRL$adjusted_p, diff_meth_TI_UCCTRL$delta_beta, "CD", "TI", "primary")
manhattan_TI_CD_organoid<-manhattan_data(rownames(diff_meth_TI_UCCTRL_organoid), diff_meth_TI_UCCTRL_organoid$p.value, diff_meth_TI_UCCTRL_organoid$adjusted_p, diff_meth_TI_UCCTRL_organoid$delta_beta, "CD", "TI", "organoid")


#In SC compare IBD vs ctrl
load(here("DNAm/data","organoid_diff_DNAm_simple_casectrl_ThreeBatch.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","organoid_diff_DNAm_simple_casectrl_ThreeBatch.RData"))
manhattan_SC_UCCD_organoid<-manhattan_data(rownames(diff_meth_SC_organoid), diff_meth_SC_organoid$p.value, diff_meth_SC_organoid$adjusted_p, diff_meth_SC_organoid$delta_beta, "IBD", "SC", "organoid")

load(here("DNAm/data","diff_DNAm_simple_casectrl.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","diff_DNAm_simple_casectrl.RData"))
manhattan_SC_UCCD_primary<-manhattan_data(rownames(diff_meth_SC), diff_meth_SC$p.value, diff_meth_SC$adjusted_p, diff_meth_SC$delta_beta, "IBD", "SC", "primary")


## EPIC only CpGs
load(file=here("DNAm/data","diff_DNAm_EPIC_only.RData"))
#load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","diff_DNAm_EPIC_only.RData"))

diff_meth_TI_UCCTRL$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL))

manhattan_TI_CD_primary_EPIC<-manhattan_data(rownames(diff_meth_TI_UCCTRL), diff_meth_TI_UCCTRL$p.value, diff_meth_TI_UCCTRL$adjusted_p, diff_meth_TI_UCCTRL$delta_beta, "CD", "TI", "primary")
manhattan_SC_UC_primary_EPIC<-manhattan_data(diff_meth_SC_grouped$CpG, diff_meth_SC_grouped$p.value_UC,diff_meth_SC_grouped$adjusted_p_UC, diff_meth_SC_grouped$db_ctrl_UC, "UC", "SC", "primary")
manhattan_SC_CD_primary_EPIC<-manhattan_data(diff_meth_SC_grouped$CpG, diff_meth_SC_grouped$p.value_CD,diff_meth_SC_grouped$adjusted_p_CD, diff_meth_SC_grouped$db_ctrl_CD, "CD", "SC", "primary")
manhattan_SC_UCCD_primary_EPIC<-manhattan_data(rownames(diff_meth_SC), diff_meth_SC$p.value,diff_meth_SC$adjusted_p, diff_meth_SC$delta_beta, "IBD", "SC", "primary")


#'# overlapping CpGs Organoid and primary
intersect(sig_hits(manhattan_SC_UCCD_primary)$CpG, sig_hits(manhattan_SC_UCCD_organoid)$CpG) #12
intersect(sig_hits(manhattan_SC_UC_primary)$CpG, sig_hits(manhattan_SC_UC_organoid)$CpG) #0
intersect(sig_hits(manhattan_SC_CD_primary)$CpG, sig_hits(manhattan_SC_CD_organoid)$CpG) #11
intersect(sig_hits(manhattan_TI_CD_primary)$CpG, sig_hits(manhattan_TI_CD_organoid)$CpG) #51

#'# overlapping CpGs primary between tests 
intersect(sig_hits(manhattan_TI_CD_primary)$CpG, sig_hits(manhattan_SC_UCCD_primary)$CpG) #0
intersect(sig_hits(manhattan_SC_UC_primary)$CpG, sig_hits(manhattan_SC_UCCD_primary)$CpG) #0
intersect(sig_hits(manhattan_SC_CD_primary)$CpG, sig_hits(manhattan_SC_UCCD_primary)$CpG) #34
intersect(sig_hits(manhattan_SC_CD_organoid)$CpG, sig_hits(manhattan_SC_UCCD_organoid)$CpG) #34
intersect(sig_hits(manhattan_SC_CD_primary)$CpG, sig_hits(manhattan_SC_UC_primary)$CpG) #0


#'# overlapping CpGs organoid between tests 
intersect(sig_hits(manhattan_TI_CD_organoid)$CpG, sig_hits(manhattan_SC_UCCD_organoid)$CpG) #837

length(intersect(sig_hits(manhattan_TI_CD_organoid)$CpG, sig_hits(manhattan_SC_CD_organoid)$CpG)) #350
length(sig_hits(manhattan_TI_CD_organoid)$CpG) #522
length(sig_hits(manhattan_SC_CD_organoid)$CpG) #7410

length(intersect(sig_hits(manhattan_TI_UC_organoid)$CpG, sig_hits(manhattan_SC_UC_organoid)$CpG)) #350
length(sig_hits(manhattan_TI_UC_organoid)$CpG) #522
length(sig_hits(manhattan_SC_UC_organoid)$CpG) #7410

intersect(sig_hits(manhattan_TI_CD_organoid)$CpG, sig_hits(manhattan_SC_UC_organoid)$CpG) #112

IBD_general<-intersect(sig_hits(manhattan_TI_CD_organoid)$CpG, sig_hits(manhattan_SC_UCCD_organoid)$CpG)
CD_specific_SC<-sig_hits(manhattan_SC_CD_organoid)$CpG
UC_specific_SC<-sig_hits(manhattan_SC_UC_organoid)$CpG

length(intersect(UC_specific_SC, CD_specific_SC)) #1025
length(intersect(intersect(UC_specific_SC, CD_specific_SC),sig_hits(manhattan_SC_UCCD_organoid)$CpG))#1023
length(intersect(UC_specific_SC, sig_hits(manhattan_SC_UCCD_organoid)$CpG)) #1302
length(UC_specific_SC)#1596
length(intersect(CD_specific_SC, sig_hits(manhattan_SC_UCCD_organoid)$CpG)) #5691
length(CD_specific_SC)#6974

length(sig_hits(manhattan_SC_UCCD_organoid)$CpG)#8277
length(sig_hits(manhattan_SC_UCCD_organoid)$CpG[which(!(sig_hits(manhattan_SC_UCCD_organoid)$CpG%in%c(CD_specific_SC, UC_specific_SC)))])#2307



length(intersect(IBD_general, CD_specific_SC)) #795

length(intersect(intersect(IBD_general, CD_specific_SC),intersect(sig_hits(manhattan_SC_UCCD_organoid)$CpG, CD_specific_SC)))#795

intersect(sig_hits(manhattan_SC_UCCD_organoid)$CpG, CD_specific_SC) #5691




#' ## Gene CpG annotation
# Made in previous scripts
EPIC_genes<-read.csv(here("DNAm/data","EPIC_ensembl_gene_annotation.csv")) # 1137194
#EPIC_genes<-read.csv(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","EPIC_ensembl_gene_annotation.csv")) # 1137194

EPIC_genes[which(EPIC_genes$IlmnID%in%intersect(sig_hits(manhattan_TI_CD_primary)$CpG, sig_hits(manhattan_TI_CD_organoid)$CpG)),]
EPIC_genes[which(EPIC_genes$IlmnID%in%intersect(sig_hits(manhattan_SC_CD_primary)$CpG, sig_hits(manhattan_SC_CD_organoid)$CpG)),]

#' ## Differentially DNAm genes
# downloaded from illumina website
anno_EPIC<-read.csv(here("DNAm/data","MethylationEPIC_v-1-0_B4.csv"), skip=7)
#anno_EPIC<-read.csv(here("/media/redgar/Seagate Portable Drive/EBI_backup/desktop_june2022/Documents/ibd/data","MethylationEPIC_v-1-0_B4.csv"), skip=7)

#' ### TI CD organoid and primary genes
TI_CD_genes<-rbind(sig_hits(manhattan_TI_CD_primary),sig_hits(manhattan_TI_CD_organoid))

TI_CD_genes$gene<-sapply(1:nrow(TI_CD_genes), function(x){
  paste0(unique(if(TI_CD_genes$CpG[x]%in%EPIC_genes$IlmnID){EPIC_genes$Gene.name[which(EPIC_genes$IlmnID==TI_CD_genes$CpG[x])]}else{"None"}), collapse = ", ")})
TI_CD_genes<-merge(TI_CD_genes, anno_EPIC[,c("IlmnID","Genome_Build","CHR","MAPINFO")], by.x="CpG", by.y="IlmnID")
TI_CD_genes<-TI_CD_genes[,c(1,7,6,5,2:4,9:11,8)]
TI_CD_genes<-TI_CD_genes[with(TI_CD_genes, order(CHR, MAPINFO)), ]
TI_CD_genes$gene<-gsub(",",";",TI_CD_genes$gene)
write.csv(TI_CD_genes, file=here("DNAm/data","CD_TI_primary_organoid_diffDNAm_gene_threebatch.csv"), row.names=F, quote = F)

TI_CD_genes_MHC1<-TI_CD_genes[grep("HLA-F|HLA-G|HLA-A|HLA-E|HLA-C|HLA-B|TAP1|TAP2|PSMB9|PSMB8|B2M|MR1|CD1D|IRF1|NLRC5", TI_CD_genes$gene),]
TI_CD_genes_MHC1$gene<-gsub(",",";",TI_CD_genes_MHC1$gene)
TI_CD_genes_MHC1
TI_CD_genes_MHC1$sample.type<-as.factor(TI_CD_genes_MHC1$sample.type)
levels(TI_CD_genes_MHC1$sample.type)<-c("organoid","IEC")
write.csv(TI_CD_genes_MHC1, file=here("DNAm/data","CD_TI_primary_organoid_diffDNAm_gene_MHC1_threebatch.csv"), row.names=F, quote = F)



#' #######
#' #'## For thesis
#' #######
#' TI_CD_genes<-rbind(sig_hits(manhattan_TI_CD_primary),sig_hits(manhattan_TI_CD_organoid))
#' 
#' TI_CD_genes$gene<-sapply(1:nrow(TI_CD_genes), function(x){
#'   paste0(unique(if(TI_CD_genes$CpG[x]%in%EPIC_genes$IlmnID){EPIC_genes$Gene.name[which(EPIC_genes$IlmnID==TI_CD_genes$CpG[x])]}else{"None"}), collapse = ", ")})
#' TI_CD_genes<-merge(TI_CD_genes, anno_EPIC[,c("IlmnID","Genome_Build","CHR","MAPINFO")], by.x="CpG", by.y="IlmnID")
#' TI_CD_genes<-TI_CD_genes[,c(1,7,6,5,2:4,9:11,8)]
#' TI_CD_genes<-TI_CD_genes[with(TI_CD_genes, order(CHR, MAPINFO)), ]
#' TI_CD_genes$gene<-gsub(",",";",TI_CD_genes$gene)
#' #write.csv(TI_CD_genes, file=here("DNAm/data","CD_TI_primary_organoid_diffDNAm_gene.csv"), row.names=F, quote = F)
#' 
#' TI_CD_genes_MHC1<-TI_CD_genes[grep("HLA-F|HLA-G|HLA-A|HLA-E|HLA-C|HLA-B|TAP1|TAP2|PSMB9|PSMB8|B2M|MR1|CD1D|IRF1|NLRC5", TI_CD_genes$gene),]
#' ## tap1 grep issue
#' TI_CD_genes_MHC1<-TI_CD_genes_MHC1[-grep("TAP12", TI_CD_genes_MHC1$gene),]
#' TI_CD_genes_MHC1$gene<-gsub(",",";",TI_CD_genes_MHC1$gene)
#' 
#' TI_CD_genes_MHC1$p.value<-scientific(TI_CD_genes_MHC1$p.value,digits=3)
#' TI_CD_genes_MHC1$fdr<-round(TI_CD_genes_MHC1$fdr,3)
#' TI_CD_genes_MHC1$db<-round(TI_CD_genes_MHC1$db,3)
#' TI_CD_genes_MHC1$sample.type<-as.factor(TI_CD_genes_MHC1$sample.type)
#' levels(TI_CD_genes_MHC1$sample.type)<-c("organoid","IEC")
#' 
#' colnames(TI_CD_genes_MHC1)<-c("CpG","Sample Type","Sample Site","Diagnosis","Delta Beta","P Value","FDR","genome_build","Chromosome","Coordinate", "Gene")
#' 
#' head(TI_CD_genes_MHC1[,c(1,9,10,11,2,6,7,5)])
#' 
#' # write.csv(TI_CD_genes_MHC1[,c(1,9,10,11,2,6,7,5)], file=here("../../thesis/redgar_PhD_thesis_main/table_data","CD_TI_primary_organoid_diffDNAm_gene_MHC1.csv"), row.names=F, quote = F)
#' # 


#' ### SC CD organoid and primary genes
SC_CD_genes<-rbind(sig_hits(manhattan_SC_CD_primary),sig_hits(manhattan_SC_CD_organoid))

SC_CD_genes$gene<-sapply(1:nrow(SC_CD_genes), function(x){
  paste0(unique(if(SC_CD_genes$CpG[x]%in%EPIC_genes$IlmnID){EPIC_genes$Gene.name[which(EPIC_genes$IlmnID==SC_CD_genes$CpG[x])]}else{"None"}), collapse = ", ")})
SC_CD_genes<-merge(SC_CD_genes, anno_EPIC[,c("IlmnID","Genome_Build","CHR","MAPINFO")], by.x="CpG", by.y="IlmnID")
SC_CD_genes<-SC_CD_genes[,c(1,7,6,5,2:4,9:11,8)]
SC_CD_genes<-SC_CD_genes[with(SC_CD_genes, order(CHR, MAPINFO)), ]
SC_CD_genes$gene<-gsub(",",";",SC_CD_genes$gene)
head(SC_CD_genes)
write.csv(SC_CD_genes, file=here("DNAm/data","CD_SC_primary_organoid_diffDNAm_gene_threebatch.csv"), row.names=F, quote = F)



#' ### SC UC organoid and primary genes
SC_UC_genes<-rbind(sig_hits(manhattan_SC_UC_primary),sig_hits(manhattan_SC_UC_organoid))

SC_UC_genes$gene<-sapply(1:nrow(SC_UC_genes), function(x){
  paste0(unique(if(SC_UC_genes$CpG[x]%in%EPIC_genes$IlmnID){EPIC_genes$Gene.name[which(EPIC_genes$IlmnID==SC_UC_genes$CpG[x])]}else{"None"}), collapse = ", ")})
SC_UC_genes<-merge(SC_UC_genes, anno_EPIC[,c("IlmnID","Genome_Build","CHR","MAPINFO")], by.x="CpG", by.y="IlmnID")
SC_UC_genes<-SC_UC_genes[,c(1,7,6,5,2:4,9:11,8)]
SC_UC_genes<-SC_UC_genes[with(SC_UC_genes, order(CHR, MAPINFO)), ]
SC_UC_genes$gene<-gsub(",",";",SC_UC_genes$gene)
head(SC_UC_genes)
write.csv(SC_UC_genes, file=here("DNAm/data","UC_SC_primary_organoid_diffDNAm_gene_threebatch.csv"), row.names=F, quote = F)


#' ### SC IBD organoid and primary genes
SC_IBD_genes<-rbind(sig_hits(manhattan_SC_UCCD_primary),sig_hits(manhattan_SC_UCCD_organoid))

SC_IBD_genes$gene<-sapply(1:nrow(SC_IBD_genes), function(x){
  paste0(unique(if(SC_IBD_genes$CpG[x]%in%EPIC_genes$IlmnID){EPIC_genes$Gene.name[which(EPIC_genes$IlmnID==SC_IBD_genes$CpG[x])]}else{"None"}), collapse = ", ")})
SC_IBD_genes<-merge(SC_IBD_genes, anno_EPIC[,c("IlmnID","Genome_Build","CHR","MAPINFO")], by.x="CpG", by.y="IlmnID")
SC_IBD_genes<-SC_IBD_genes[,c(1,7,6,5,2:4,9:11,8)]
SC_IBD_genes<-SC_IBD_genes[with(SC_IBD_genes, order(CHR, MAPINFO)), ]
SC_IBD_genes$gene<-gsub(",",";",SC_IBD_genes$gene)
head(SC_IBD_genes)
write.csv(SC_IBD_genes, file=here("DNAm/data","IBD_SC_primary_organoid_diffDNAm_gene_threebatch.csv"), row.names=F, quote = F)


#' ### SC CD organoid and primary genes
SC_CD_genes<-rbind(sig_hits(manhattan_SC_CD_primary),sig_hits(manhattan_SC_CD_organoid))

SC_CD_genes$gene<-sapply(1:nrow(SC_CD_genes), function(x){
  paste0(unique(if(SC_CD_genes$CpG[x]%in%EPIC_genes$IlmnID){EPIC_genes$Gene.name[which(EPIC_genes$IlmnID==SC_CD_genes$CpG[x])]}else{"None"}), collapse = ", ")})
SC_CD_genes<-merge(SC_CD_genes, anno_EPIC[,c("IlmnID","Genome_Build","CHR","MAPINFO")], by.x="CpG", by.y="IlmnID")
SC_CD_genes<-SC_CD_genes[,c(1,7,6,5,2:4,9:11,8)]
SC_CD_genes<-SC_CD_genes[with(SC_CD_genes, order(CHR, MAPINFO)), ]
SC_CD_genes$gene<-gsub(",",";",SC_CD_genes$gene)
head(SC_CD_genes)
write.csv(SC_CD_genes, file=here("DNAm/data","CD_SC_primary_organoid_diffDNAm_gene_threebatch.csv"), row.names=F, quote = F)


#'# Gene level overlap
#' # IBD genes
sig_CpG<-rbind(sig_hits(manhattan_TI_UC_primary), sig_hits(manhattan_SC_UC_primary), sig_hits(manhattan_SC_CD_primary), 
               sig_hits(manhattan_TI_UC_organoid), sig_hits(manhattan_SC_UC_organoid), sig_hits(manhattan_SC_CD_organoid), 
               sig_hits(manhattan_TI_CD_primary),sig_hits(manhattan_TI_CD_organoid),
               sig_hits(manhattan_TI_CD_primary_EPIC),sig_hits(manhattan_SC_UC_primary_EPIC),sig_hits(manhattan_SC_UCCD_primary_EPIC),
               sig_hits(manhattan_SC_UCCD_primary),sig_hits(manhattan_SC_UCCD_organoid))


IBD_genes<-merge(sig_CpG, EPIC_genes[,c(3,6,7,9)], by.x="CpG", by.y="IlmnID")

tapply(sig_CpG$CpG, list(sig_CpG$tissue, sig_CpG$diagnosis, sig_CpG$sample.type), function(x) length(x))

CD_TI<-IBD_genes[which(IBD_genes$tissue=="TI" & IBD_genes$diagnosis=="CD"),]
length(unique(CD_TI$Gene.name))
CD_SC<-IBD_genes[which(IBD_genes$tissue=="SC" & IBD_genes$diagnosis=="CD"),]
length(unique(CD_SC$Gene.name))
UC_SC<-IBD_genes[which(IBD_genes$tissue=="SC" & IBD_genes$diagnosis=="UC"),]
length(unique(UC_SC$Gene.name))
IBD_SC<-IBD_genes[which(IBD_genes$tissue=="SC" & IBD_genes$diagnosis=="IBD"),]
length(unique(IBD_SC$Gene.name))

CD_TI_primary<-unique(CD_TI$Gene.name[which(CD_TI$sample.type=="primary")])
CD_TI_organoid<-unique(CD_TI$Gene.name[which(CD_TI$sample.type=="organoid")])

CD_SC_primary<-unique(CD_SC$Gene.name[which(CD_SC$sample.type=="primary")])
CD_SC_organoid<-unique(CD_SC$Gene.name[which(CD_SC$sample.type=="organoid")])

UC_SC_primary<-unique(UC_SC$Gene.name[which(UC_SC$sample.type=="primary")])
UC_SC_organoid<-unique(UC_SC$Gene.name[which(UC_SC$sample.type=="organoid")])

IBD_SC_primary<-unique(IBD_SC$Gene.name[which(IBD_SC$sample.type=="primary")])
IBD_SC_organoid<-unique(IBD_SC$Gene.name[which(IBD_SC$sample.type=="organoid")])


length(intersect(CD_TI_primary, IBD_SC_primary))#0
length(intersect(CD_TI_organoid, IBD_SC_organoid))#216

length(intersect(CD_SC_primary, IBD_SC_primary))#31
length(intersect(CD_SC_organoid, IBD_SC_organoid))#548

length(intersect(UC_SC_primary, IBD_SC_primary))#2
length(intersect(UC_SC_organoid, IBD_SC_organoid))#49

length(intersect(CD_TI_organoid, CD_SC_organoid))#340
length(CD_TI_organoid)#422
length(CD_SC_organoid)#4324

TI_CD_genes_MHC1<-CD_TI_organoid[grep("HLA-F|HLA-G|HLA-A|HLA-E|HLA-C|HLA-B|TAP1|TAP2|PSMB9|PSMB8|B2M|IRF1|NLRC5", CD_TI_organoid)]
CD_SC_organoid_MHCI<-CD_SC_organoid[grep("HLA-F|HLA-G|HLA-A|HLA-E|HLA-C|HLA-B|TAP1|TAP2|PSMB9|PSMB8|B2M|IRF1|NLRC5", CD_SC_organoid)]



#' MHCI
HLA_IBD<-IBD_genes[grep("HLA",IBD_genes$Gene.name),]
tapply(HLA_IBD$Gene.name, list(HLA_IBD$tissue, HLA_IBD$diagnosis, HLA_IBD$sample.type), function(x) length(unique(x)))
HLA_IBD[which(HLA_IBD$sample.type=="primary"),]
# HLA-E in both

MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","MR1","CD1D","IRF1","NLRC5")
MHCI_diff<-IBD_genes[which(IBD_genes$Gene.name%in%c(MHCI,"TAPBP","CIITA")),]
unique(MHCI_diff[which(MHCI_diff$sample.type=="primary"),"Gene.name"])
unique(MHCI_diff[which(MHCI_diff$sample.type=="organoid"),"Gene.name"])

IBD_SC[which(IBD_SC$Gene.name%in%c(MHCI,"TAPBP","CIITA")),]

#' ## delta beta plot 
db_plot<-function(manhattanprimary,manhattanorganoid, tissue, diagnosis){
  plot_df<-merge(manhattanprimary,manhattanorganoid, by="CpG")
  print(cor(plot_df$db.x, plot_df$db.y, method="spearman"))
  
  plot_df$color<-sapply(1:nrow(plot_df), function(x){
    if(is.na(plot_df$fdr.x[x])){NA}else{
    if(plot_df$fdr.x[x] < 0.05 & abs(plot_df$db.x[x])>0.05 & plot_df$fdr.y[x] < 0.05 & abs(plot_df$db.y[x])>0.05){"Primary and\nOrganoid"}else{
      if(plot_df$fdr.x[x] < 0.05 & abs(plot_df$db.x[x])>0.05){"Primary"}else{
        if(plot_df$fdr.y[x] < 0.05 & abs(plot_df$db.y[x])>0.05){"Organoid"}else{"Not\nSignificant"}
      }
    }
  }})
  
  plot_df<-plot_df[(order(plot_df$color)),]
  print(paste("at the",nrow(plot_df[which(plot_df$color=="Primary and\nOrganoid"),]), "CpGs sig in both the cor is ", cor(plot_df$db.x[which(plot_df$color=="Primary and\nOrganoid")], plot_df$db.y[which(plot_df$color=="Primary and\nOrganoid")], method="spearman")))
  
  ggplot(plot_df, aes(db.x,db.y))+
    geom_point(aes(color=color, alpha=color), shape=19,size=2)+theme_bw()+th+
    geom_hline(yintercept=0, color="grey30")+
    geom_vline(xintercept=0, color="grey30")+
    geom_vline(xintercept=c(0.05, -0.05), color="grey70")+
    geom_hline(yintercept=c(0.05, -0.05), color="grey70")+
    xlim(-0.35, 0.35)+ylim(-0.35, 0.35)+
    #scale_color_manual(values=c("lightgrey","#01665e","#b2182b","#7a0177"), name="Significant")+#1874CD
    scale_color_manual(values=c("lightgrey","#65b3ff","#99ccff","#1874CD"), name="Significant")+
    scale_alpha_manual(values=c(0.25,1,1,1), name="Significant")+
    stat_smooth(method="lm", color="black")+
    xlab(paste(tissue, diagnosis, "Primary Delta Beta"))+ ylab(paste(tissue, diagnosis, "Organoid Delta Beta"))+
    th
}

db_plot(manhattan_SC_UC_primary,manhattan_SC_UC_organoid,"SC","UC")
ggsave(here("DNAm/figs/jpeg","delta_beta_cor_SCUC_threebatch.jpeg"),  width = 9, height = 7.5)
# ggsave(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/figs/jpeg","delta_beta_cor_SCUC_threebatch.jpeg"),  width = 9, height = 7.5)
# ggsave(here("/home/redgar/Documents/EBI/MHCI/IBD_bulk_integration/DNAm/figs/jpeg","delta_beta_cor_SCUC_threebatch.jpeg"),  width = 9, height = 7.5)

db_plot(manhattan_SC_CD_primary,manhattan_SC_CD_organoid,"SC","CD") 
ggsave(here("DNAm/figs/jpeg","delta_beta_cor_SCCD_threebatch.jpeg"),  width = 9, height = 7.5)
# ggsave(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/figs/jpeg","delta_beta_cor_SCCD_threebatch.jpeg"),  width = 9, height = 7.5)
# ggsave(here("/home/redgar/Documents/EBI/MHCI/IBD_bulk_integration/DNAm/figs/jpeg","delta_beta_cor_SCCD_threebatch.jpeg"),  width = 9, height = 7.5)

db_plot(manhattan_TI_CD_primary,manhattan_TI_CD_organoid,"TI","CD") 
ggsave(here("DNAm/figs/jpeg","delta_beta_cor_TICD_threebatch.jpeg"),  width = 6, height = 4.5)
# ggsave(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/figs/jpeg","delta_beta_cor_TICD_threebatch.jpeg"),  width = 6, height = 4.5)
# ggsave(here("/home/redgar/Documents/EBI/MHCI/IBD_bulk_integration/DNAm/figs/jpeg","delta_beta_cor_TICD_threebatch.jpeg"),  width = 9, height = 7.5)

db_plot(manhattan_SC_UCCD_primary,manhattan_SC_UCCD_organoid,"SC","IBD") 
ggsave(here("DNAm/figs/jpeg","delta_beta_cor_SCIBD_threebatch.jpeg"),  width = 6, height = 4.5)
# ggsave(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/figs/jpeg","delta_beta_cor_SCIBD_threebatch.jpeg"),  width = 6, height = 4.5)
# ggsave(here("/home/redgar/Documents/EBI/MHCI/IBD_bulk_integration/DNAm/figs/jpeg","delta_beta_cor_SCIBD_threebatch.jpeg"),  width = 9, height = 7.5)

#' #######
#' #'## For thesis
#' #######
#' db_plot(manhattan_TI_CD_primary,manhattan_TI_CD_organoid,"TI","CD") 
#' ggsave(here("../../thesis/redgar_PhD_thesis_main/Chapter3/Figs/Raster","delta_beta_cor_TICD.jpeg"),  width = 6, height = 4.5)
#' 
#' #write.csv(TI_CD_genes_MHC1[,c(1,9,10,11,2,6,7,5)], file=here("../../thesis/redgar_PhD_thesis_main/table_data","CD_TI_primary_organoid_diffDNAm_gene_MHC1.csv"), row.names=F, quote = F)
#' 
#' 


#' ## for MUC1 focused presentation
CpG_OI<-c("cg24512973","cg18804777")

plot_df<-merge(manhattan_TI_CD_primary,manhattan_TI_CD_organoid, by="CpG")
print(cor(plot_df$db.x, plot_df$db.y, method="spearman"))

plot_df$color<-sapply(1:nrow(plot_df), function(x){
  if(plot_df$fdr.x[x] < 0.05 & abs(plot_df$db.x[x])>0.05 & plot_df$fdr.y[x] < 0.05 & abs(plot_df$db.y[x])>0.05){"Primary and\nOrganoid"}else{
    if(plot_df$fdr.x[x] < 0.05 & abs(plot_df$db.x[x])>0.05){"Primary"}else{
      if(plot_df$fdr.y[x] < 0.05 & abs(plot_df$db.y[x])>0.05){"Organoid"}else{"Not\nSignificant"}
    }
  }
})

plot_df<-plot_df[(order(plot_df$color)),]

plot_df$hightlight<-" "
plot_df$hightlight[which(plot_df$CpG%in%CpG_OI)]<-"MUC1 CpG"

ggplot(plot_df, aes(db.x,db.y))+
  geom_point(aes(fill=color, alpha=color, color=hightlight), shape=21,size=2)+theme_bw()+th+
  geom_hline(yintercept=0, color="grey30")+
  geom_vline(xintercept=0, color="grey30")+
  geom_vline(xintercept=c(0.05, -0.05), color="grey70")+
  geom_hline(yintercept=c(0.05, -0.05), color="grey70")+
  xlim(-0.3, 0.3)+ylim(-0.3, 0.3)+
  scale_fill_manual(values=c("lightgrey","#65b3ff","#99ccff","#1874CD"), name="Differential\nDNAm")+
  scale_alpha_manual(values=c(0.25,1,1,1), name="Differential\nDNAm")+
  scale_color_manual(values=c("white","black"), name="")+
  stat_smooth(method="lm", color="black")+
  xlab("Primary Delta Beta")+ ylab("Organoid Delta Beta")+th

ggsave(here("DNAm/figs/jpeg","delta_beta_cor_TICD_MUC1_threebatch.jpeg"),  width = 6, height = 4.5)







###############
#' # Function to plot any gene or CpG set
###############
#' ### Organoid
load(file=paste(here("DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))
# load(file=paste(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))

## fix 202 and 203
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605870061_R03C01")]<-"hold"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605880096_R04C01")]<-"205605870061_R03C01"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="hold")]<-"205605880096_R04C01"



#'### Primary
#' remove rescopes 
sampleinfo_og<-sampleinfo_DNAm[which(is.na(sampleinfo_DNAm$passage.or.rescope.no)),]
ibd_combo_og<-combat_ibd_Beta[,which(colnames(combat_ibd_Beta)%in%sampleinfo_og$array.id)]
identical(colnames(ibd_combo_og),sampleinfo_og$array.id)
sampleinfo_og$diagnosis<-as.factor(sampleinfo_og$diagnosis_grouped)

#' Just EPIC
load(here("DNAm/data/ibd_adjusted_betas_EPIC_only.RData"))
#load(here("/media/redgar/Seagate Portable Drive/EBI_backup/desktop_august2022/Documents/ibd/IBD_bulk_integration/DNAm/data/ibd_adjusted_betas_EPIC_only.RData"))

sampleinfo_DNAm_EPIC<-sampleinfo_DNAm[which(sampleinfo_DNAm$array.type=="EPIC"),]
identical(sampleinfo_DNAm_EPIC$array.id, colnames(ibd_epic_adjusted))
## only need for sites not tested on both arrays
ibd_beta_epic<-ibd_epic_adjusted[which(!(rownames(ibd_epic_adjusted)%in%rownames(ibd_combo_og))),]


######################
#' ## Gene plot organoid and primary
######################
all_stats<-rbind(manhattan_TI_UC_primary, manhattan_SC_UC_primary, manhattan_SC_CD_primary, 
                 manhattan_TI_UC_organoid, manhattan_SC_UC_organoid, manhattan_SC_CD_organoid, 
                 manhattan_TI_CD_primary,manhattan_TI_CD_organoid,
                 manhattan_TI_CD_primary_EPIC,manhattan_SC_UC_primary_EPIC,manhattan_SC_UCCD_primary_EPIC,
                 manhattan_SC_UCCD_primary, manhattan_SC_UCCD_organoid)

meta_merge<-rbind(sampleinfo_og[,c("array.id","sample.site","diagnosis")], epic.organoid_combined[,c("array.id","sample.site","diagnosis")])


DMR_function<-function(input, gene_name,around_sig,plot_type){#, plot_type
  if(input=="CpG"){ CpG_OI<-gene_name}
  if(input=="gene"){ CpG_OI<-EPIC_genes[which(EPIC_genes$Gene.name==gene_name),]$IlmnID}
  
  all_stats_OI<-all_stats[which(all_stats$CpG%in%CpG_OI),]
  sig<-all_stats_OI[which(all_stats_OI$fdr<0.05 & abs(all_stats_OI$db)>=0.05),]
  sig<-merge(sig, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
  print(sig)
  
  #CpG_OI<-CpG_OI[which(CpG_OI%in%sig$CpG)]
  if(nrow(sig)==0){print("No Sig. IBD CpG")}else{
    
    if(!is.na(around_sig) & input=="gene"){
      CpG_OI<-anno_EPIC[which(anno_EPIC$MAPINFO>=min(sig$MAPINFO)-around_sig & anno_EPIC$MAPINFO<=max(sig$MAPINFO)+around_sig & anno_EPIC$CHR==unique(EPIC_genes[which(EPIC_genes$Gene.name==gene_name),"CHR"])),"IlmnID"]
      all_stats_OI<-all_stats[which(all_stats$CpG%in%CpG_OI),]
      sig<-all_stats_OI[which(all_stats_OI$fdr <0.05 & abs(all_stats_OI$db)>=0.05),]
      sig<-merge(sig, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
    }}
  
  
  
  betas_org<-as.data.frame(combat_organoid_Beta[which(rownames(combat_organoid_Beta)%in%CpG_OI),])
  betas_org$CpG<-rownames(betas_org)
  betas_org<-melt(betas_org)
  betas_org$sample.type<-"organoid"
  betas_org$sample.n<-"full"
  
  betas_pri<-as.data.frame(ibd_combo_og[which(rownames(ibd_combo_og)%in%CpG_OI),])
  if(ncol(betas_pri)==1){
    betas_pri$variable<-rownames(betas_pri)
    betas_pri$CpG<-rownames(ibd_combo_og)[which(rownames(ibd_combo_og)%in%CpG_OI)]
    colnames(betas_pri)<-c("value","variable","CpG")
    betas_pri<-betas_pri[,c("CpG","variable","value")]
  }else{
    betas_pri$CpG<-rownames(betas_pri)
    betas_pri<-melt(betas_pri)}
  betas_pri$sample.type<-"primary"
  betas_pri$sample.n<-"full"
  
  betas_pri_epic<-as.data.frame(ibd_beta_epic[which(rownames(ibd_beta_epic)%in%CpG_OI),])
  
  if(length(which(rownames(ibd_beta_epic)%in%CpG_OI))==0){betas_pri_epic[1,]<-NA}else{
    if(length(which(rownames(ibd_beta_epic)%in%CpG_OI))==1){
      colnames(betas_pri_epic)<-rownames(ibd_beta_epic)[which(rownames(ibd_beta_epic)%in%CpG_OI)]
      betas_pri_epic<-as.data.frame(t(betas_pri_epic))}}
  
  betas_pri_epic$CpG<-rownames(betas_pri_epic)
  betas_pri_epic<-melt(betas_pri_epic)
  betas_pri_epic$sample.type<-"primary"
  betas_pri_epic$sample.n<-"sub"
  
  betas<-rbind(betas_org,betas_pri,betas_pri_epic)
  
  betas<-merge(betas, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
  plt<-merge(betas, meta_merge, by.x="variable", by.y="array.id")
  
  plt<-plt[which(plt$sample.site!="AC"),]
  plt<-plt[which(plt$sample.site!="Duo"),]
  
  
  sum_stat<-melt(tapply(plt$value, list(plt$MAPINFO,plt$diagnosis,plt$sample.site, plt$sample.type), mean))
  colnames(sum_stat)<-c("MAPINFO","diagnosis","sample.site","sample.type","value")
  sum_stat<-sum_stat[which(!(is.na(sum_stat$value))),]
  
  colnames(sig)[6]<-"sample.site"
  
  
  if(plot_type=="dmr"){
    print(ggplot()+ #geom_rect(aes(xmin = MAPINFO-stat_highlight, xmax = MAPINFO+stat_highlight, ymin = 0, ymax = 1), fill = '#f2f500', alpha = 1, data = sig)+
            geom_vline(data = sig,
                       aes(xintercept = MAPINFO),color = 'black', size=1)+
            geom_point(aes(MAPINFO, value, fill=diagnosis, alpha=as.factor(sample.n)), plt, size=1.5,shape=21, color="white")+
            geom_line(aes(MAPINFO,value, color=diagnosis),sum_stat,size=0.5)+
            facet_grid(sample.type~sample.site)+
            colscale_diagnosis+fillscale_diagnosis+theme_bw()+th+ylim(0,1)+ylab("DNAm Beta Value")+
            scale_alpha_manual(values=c(1,0.15),guide=F))}else{
              
              if(plot_type=="box"){
                plt$new_color<-paste(plt$diagnosis, plt$sample.n)
                plt$new_color<-factor(plt$new_color, levels=c("CD full","CD sub","UC full","UC sub","IBD-U full","Control full", "Control sub"))
                
                
                p<-ggplot()+geom_boxplot(aes(CpG, value,fill=new_color, color=sample.n), plt, outlier.size = 0.75)+facet_grid(sample.type~sample.site)+
                  scale_color_manual(values=c("black","grey"))+fillscale_diagnosis_sample+theme_bw()+th+ylim(0,1)+ylab("DNAm Beta Value")+
                  theme(plot.margin = margin(0.5, 0.1, 0.1, 0.5, "cm"),
                        axis.text = element_text(size =10, angle=45,color="black", hjust=1),
                        axis.title = element_text(size =12),
                        legend.text = element_text(size =12),
                        legend.title = element_text(size =12),
                        strip.text.x = element_text(size = 12))
                
                if(nrow(sig)>0){
                  # plt<-plt[which(plt$CpG%in%sig$CpG),]
                  # sig_plt<-plt[,c(2,4,7,8)][!duplicated(plt[,c(2,4,7,8)]),]
                  sig$sig<-"*"
                  # sig_plt<-merge(sig_plt, sig[,c(1,5:7,9)],by=c("CpG","sample.type", "sample.site", "diagnosis"),all.x=T)
                  # p+geom_text(aes(CpG,label = sig, y=0.95),data=sig_plt, size=8)}else{p}}
                  p+geom_text(aes(CpG,label = paste(sig, diagnosis), y=0.95),data=sig, size=3)}else{p}}
            }
}

            # input="CpG"
            # gene_name=c("cg07839457","cg07862320","cg16411857","cg18440314")
            # around_sig=0
            # plot_type="box"
#' # Plot Key genes
DMR_function("gene","NLRC5", 0, "box")
DMR_function("CpG",c("cg07839457","cg07862320","cg16411857","cg18440314"), 0, "box")
DMR_function("gene","NLRC5", 0, "dmr")

pltsave=DMR_function("gene","NLRC5", 0, "dmr")
pltsave
ggsave(pltsave, file=here("DNAm/figs","NLRC5_DNAm_threebatch.pdf"), width=8, height=5)


pltsave=DMR_function("CpG",c("cg07839457","cg07862320"), 0, "box")
pltsave
ggsave(pltsave, file=here("DNAm/figs","NLRC5_DNAm_sig_threebatch.pdf"), width=8, height=5)
ggsave(pltsave, file=here("DNAm/figs/jpeg","NLRC5_DNAm_sig_threebatch.jpeg"), width=8, height=5)


pltsave=DMR_function("CpG",c("cg11706729","cg02756056"), 0, "box")
pltsave
ggsave(pltsave, file=here("DNAm/figs","TAP1_DNAm_threebatch.pdf"), width=8, height=5)
ggsave(pltsave, file=here("DNAm/figs/jpeg","TAP1_DNAm_threebatch.jpeg"), width=8, height=5)


pltsave=DMR_function("gene","PSMB8", 0, "dmr")
pltsave
ggsave(pltsave, file=here("DNAm/figs","PSMB8_DNAm_threebatch.pdf"), width=8, height=5)

pltsave=DMR_function("gene","PLA2G2A", 5000, "dmr")
pltsave
ggsave(pltsave, file=here("DNAm/figs","PLA2G2A_DNAm_threebatch.pdf"), width=8, height=5)

pltsave=DMR_function("gene","B2M", 5000, "dmr")
pltsave
ggsave(pltsave, file=here("DNAm/figs","B2M_DNAm_threebatch.pdf"), width=8, height=5)


pltsave=DMR_function("CpG",c("cg24512973","cg18804777"), 0, "box")
pltsave
ggsave(pltsave, file=here("DNAm/figs","MUC1_DNAm_threebatch.pdf"), width=8, height=5)


pltsave=DMR_function("CpG",c("cg03783907", "cg12537337", "cg12882188", "cg15460604"), 0, "box")
pltsave
ggsave(pltsave, file=here("DNAm/figs","BTN_DNAm_sig_threebatch.pdf"), width=8, height=5)
ggsave(pltsave, file=here("DNAm/figs/jpeg","BTN_DNAm_sig_threebatch.jpeg"), width=8, height=5)

print("START OF MAIN BOXPLOTS")


#' # NLRC5 to present
CpG_OI<-c("cg07839457","cg07862320")

all_stats_OI<-all_stats[which(all_stats$CpG%in%CpG_OI),]
sig<-all_stats_OI[which(all_stats_OI$fdr<0.05 & abs(all_stats_OI$db)>=0.05),]
sig<-merge(sig, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")

betas_org<-as.data.frame(combat_organoid_Beta[which(rownames(combat_organoid_Beta)%in%CpG_OI),])
betas_org$CpG<-rownames(betas_org)
betas_org<-melt(betas_org)
betas_org$sample.type<-"organoid"
betas_org$sample.n<-"full"

betas_pri<-as.data.frame(ibd_combo_og[which(rownames(ibd_combo_og)%in%CpG_OI),])
betas_pri$variable<-rownames(betas_pri)
betas_pri$CpG<-rownames(ibd_combo_og)[which(rownames(ibd_combo_og)%in%CpG_OI)]
colnames(betas_pri)<-c("value","variable","CpG")
betas_pri<-betas_pri[,c("CpG","variable","value")]
betas_pri$sample.type<-"primary"
betas_pri$sample.n<-"full"

betas_pri_epic<-as.data.frame(ibd_beta_epic[which(rownames(ibd_beta_epic)%in%CpG_OI),])
if(length(which(rownames(ibd_beta_epic)%in%CpG_OI))==1){
  colnames(betas_pri_epic)<-rownames(ibd_beta_epic)[which(rownames(ibd_beta_epic)%in%CpG_OI)]
  betas_pri_epic<-as.data.frame(t(betas_pri_epic))
}
betas_pri_epic$CpG<-rownames(betas_pri_epic)
betas_pri_epic<-melt(betas_pri_epic)
betas_pri_epic$sample.type<-"primary"
betas_pri_epic$sample.n<-"sub"

betas<-rbind(betas_org,betas_pri,betas_pri_epic)

betas<-merge(betas, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
plt<-merge(betas, meta_merge, by.x="variable", by.y="array.id")

plt<-plt[which(plt$sample.site=="TI"),]

sum_stat<-melt(tapply(plt$value, list(plt$MAPINFO,plt$diagnosis,plt$sample.site, plt$sample.type), mean))
colnames(sum_stat)<-c("MAPINFO","diagnosis","sample.site","sample.type","value")
sum_stat<-sum_stat[which(!(is.na(sum_stat$value))),]

colnames(sig)[6]<-"sample.site"


plt$new_color<-paste(plt$diagnosis, plt$sample.n)
plt<-plt[which(plt$CpG%in%sig$CpG),]

###rm IBD-U from TI
plt<-plt[which(!(plt$variable%in%epic.organoid_exclude_other$array.id[which(epic.organoid_exclude_other$diagnosis=="IBD-U")])),]

plt$diagnosis<-factor(plt$diagnosis, levels=c("Control","UC","CD"))

p<-ggplot()+geom_boxplot(aes(diagnosis, value,fill=new_color, color=sample.n), plt, outlier.shape = NA)+facet_grid(sample.type~CpG)+
  geom_jitter(aes(diagnosis, value,fill=new_color, color=sample.n), plt, width=0.25, shape=21, size=2)+
  scale_color_manual(values=c("black","grey40"))+fillscale_diagnosis_sample+theme_bw()+th+ylim(0,1)+ylab("DNAm Beta Value")+xlab("Diagnosis")+
  theme(legend.position="none",
        axis.title = element_text(size =12),
        legend.text = element_text(size =12),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))

p
sig_plt<-plt[,c(2,4,7,8)][!duplicated(plt[,c(2,4,7,8)]),]
sig$sig<-"*"
print(sig)
sig_plt<-merge(sig_plt, sig[,c(1,5:7,9)],by=c("CpG","sample.type","diagnosis"),all.x=T)
p+geom_text(aes(diagnosis,label = sig, y=0.95),data=sig_plt, size=8)

ggsave(file=here("DNAm/figs","NLRC5_DNAm_box_threebatch.pdf"), width=7, height=5)
ggsave(file=here("DNAm/figs/jpeg","NLRC5_DNAm_box_threebatch.jpeg"), width=7, height=5)



CpG_OI<-c("cg07839457","cg07862320")

all_stats_OI<-all_stats[which(all_stats$CpG%in%CpG_OI),]
sig<-all_stats_OI[which(all_stats_OI$fdr<0.05 & abs(all_stats_OI$db)>=0.05),]
sig<-merge(sig, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")

betas_org<-as.data.frame(combat_organoid_Beta[which(rownames(combat_organoid_Beta)%in%CpG_OI),])
betas_org$CpG<-rownames(betas_org)
betas_org<-melt(betas_org)
betas_org$sample.type<-"organoid"
betas_org$sample.n<-"full"

betas_pri<-as.data.frame(ibd_combo_og[which(rownames(ibd_combo_og)%in%CpG_OI),])
betas_pri$variable<-rownames(betas_pri)
betas_pri$CpG<-rownames(ibd_combo_og)[which(rownames(ibd_combo_og)%in%CpG_OI)]
colnames(betas_pri)<-c("value","variable","CpG")
betas_pri<-betas_pri[,c("CpG","variable","value")]
betas_pri$sample.type<-"primary"
betas_pri$sample.n<-"full"

betas_pri_epic<-as.data.frame(ibd_beta_epic[which(rownames(ibd_beta_epic)%in%CpG_OI),])
if(length(which(rownames(ibd_beta_epic)%in%CpG_OI))==1){
  colnames(betas_pri_epic)<-rownames(ibd_beta_epic)[which(rownames(ibd_beta_epic)%in%CpG_OI)]
  betas_pri_epic<-as.data.frame(t(betas_pri_epic))
}
betas_pri_epic$CpG<-rownames(betas_pri_epic)
betas_pri_epic<-melt(betas_pri_epic)
betas_pri_epic$sample.type<-"primary"
betas_pri_epic$sample.n<-"sub"

betas<-rbind(betas_org,betas_pri,betas_pri_epic)

betas<-merge(betas, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
plt<-merge(betas, meta_merge, by.x="variable", by.y="array.id")

plt<-plt[which(plt$sample.site=="SC"),]

sum_stat<-melt(tapply(plt$value, list(plt$MAPINFO,plt$diagnosis,plt$sample.site, plt$sample.type), mean))
colnames(sum_stat)<-c("MAPINFO","diagnosis","sample.site","sample.type","value")
sum_stat<-sum_stat[which(!(is.na(sum_stat$value))),]

colnames(sig)[6]<-"sample.site"


plt$new_color<-paste(plt$diagnosis, plt$sample.n)
plt<-plt[which(plt$CpG%in%sig$CpG),]

plt$diagnosis<-factor(plt$diagnosis, levels=c("Control","IBD-U","UC","CD"))

p<-ggplot()+geom_boxplot(aes(diagnosis, value,fill=new_color, color=sample.n), plt, outlier.shape = NA)+facet_grid(sample.type~CpG)+
  geom_jitter(aes(diagnosis, value,fill=new_color, color=sample.n), plt, width=0.25, shape=21, size=2)+
  scale_color_manual(values=c("black","grey40"))+fillscale_diagnosis_sample+theme_bw()+th+ylim(0,1)+ylab("DNAm Beta Value")+xlab("Diagnosis")+
  theme(legend.position="none",plot.margin = margin(0.5, 0.1, 0.1, 0.5, "cm"),
        axis.text = element_text(size =10, angle=45,color="black", hjust=1),
        axis.title = element_text(size =12),
        legend.text = element_text(size =12),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))

p
sig<-sig[which(sig$sample.site=="SC"),]
sig_plt<-plt[,c(2,4,7,8)][!duplicated(plt[,c(2,4,7,8)]),]
sig$sig<-"*"
sig$diagnosis[which(sig$diagnosis=="IBD")]<-"Control"
sig_plt<-merge(sig_plt, sig[,c(1,5:7,9)],by=c("CpG","sample.type","diagnosis"),all.x=T)
p+geom_text(aes(diagnosis,label = sig, y=0.95),data=sig_plt, size=8)

ggsave(file=here("DNAm/figs","NLRC5_DNAm_box_SC_threebatch.pdf"),  width=7, height=5)
ggsave(file=here("DNAm/figs/jpeg","NLRC5_DNAm_box_SC_threebatch.jpeg"), width=7, height=5)



#' # TAP1 to present

gene_name="TAP1"
sample.site="TI"

CpG_OI<-c("cg11706729","cg02756056")

all_stats_OI<-all_stats[which(all_stats$CpG%in%CpG_OI),]
sig<-all_stats_OI[which(all_stats_OI$fdr<0.05 & abs(all_stats_OI$db)>=0.05),]
sig<-merge(sig, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")

betas_org<-as.data.frame(combat_organoid_Beta[which(rownames(combat_organoid_Beta)%in%CpG_OI),])
betas_org$CpG<-rownames(betas_org)
betas_org<-melt(betas_org)
betas_org$sample.type<-"organoid"
betas_org$sample.n<-"full"

betas_pri<-as.data.frame(ibd_combo_og[which(rownames(ibd_combo_og)%in%CpG_OI),])
betas_pri$CpG<-rownames(betas_pri)
betas_pri<-melt(betas_pri)
betas_pri$sample.type<-"primary"
betas_pri$sample.n<-"full"

betas_pri_epic<-as.data.frame(ibd_beta_epic[which(rownames(ibd_beta_epic)%in%CpG_OI),])
if(length(which(rownames(ibd_beta_epic)%in%CpG_OI))==1){
  colnames(betas_pri_epic)<-rownames(ibd_beta_epic)[which(rownames(ibd_beta_epic)%in%CpG_OI)]
  betas_pri_epic<-as.data.frame(t(betas_pri_epic))
  betas_pri_epic$CpG<-rownames(betas_pri_epic)
  betas_pri_epic<-melt(betas_pri_epic)
  betas_pri_epic$sample.type<-"primary"
  betas_pri_epic$sample.n<-"sub"
  
  betas<-rbind(betas_org,betas_pri,betas_pri_epic)}else{betas<-rbind(betas_org,betas_pri)}



betas<-merge(betas, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
plt<-merge(betas, meta_merge, by.x="variable", by.y="array.id")

plt<-plt[which(plt$sample.site=="TI"),]
###rm IBD-U from TI
plt<-plt[which(!(plt$variable%in%epic.organoid_exclude_other$array.id[which(epic.organoid_exclude_other$diagnosis=="IBD-U")])),]

sum_stat<-melt(tapply(plt$value, list(plt$MAPINFO,plt$diagnosis,plt$sample.site, plt$sample.type), mean))
colnames(sum_stat)<-c("MAPINFO","diagnosis","sample.site","sample.type","value")
sum_stat<-sum_stat[which(!(is.na(sum_stat$value))),]

colnames(sig)[6]<-"sample.site"


plt$new_color<-paste(plt$diagnosis, plt$sample.n)
plt<-plt[which(plt$CpG%in%sig$CpG),]

plt$diagnosis<-factor(plt$diagnosis, levels=c("Control","UC","CD"))

p<-ggplot()+geom_boxplot(aes(diagnosis, value,fill=new_color, color=sample.n), plt, outlier.shape = NA)+facet_grid(sample.type~CpG)+
  geom_jitter(aes(diagnosis, value,fill=new_color, color=sample.n), plt, width=0.25, shape=21, size=2)+
  scale_color_manual(values=c("black","grey40"))+fillscale_diagnosis_sample+theme_bw()+th+ylim(0,1)+ylab("DNAm Beta Value")+xlab("Diagnosis")+
  theme(legend.position="none",
        axis.text = element_text(size =12, color="black"),
        axis.title = element_text(size =12),
        legend.text = element_text(size =12),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))

p
sig_plt<-plt[,c(2,4,7,8)][!duplicated(plt[,c(2,4,7,8)]),]
sig$sig<-"*"
sig_plt<-merge(sig_plt, sig[,c(1,5:7,9)],by=c("CpG","sample.type","diagnosis"),all.x=T)
p+geom_text(aes(diagnosis,label = sig, y=0.95),data=sig_plt, size=8)

ggsave(file=here("DNAm/figs","TAP1_TI_DNAm_box_threebatch.pdf"), width=7, height=5)
ggsave(file=here("DNAm/figs/jpeg","TAP1_TI_DNAm_box_threebatch.jpeg"), width=7, height=5)

print("BOX PLOTS UP TO DATE")




#'#### Celiac comparison figure
CpG_OI<-c("cg07839457","cg26033526")

all_stats_OI<-all_stats[which(all_stats$CpG%in%CpG_OI),]
sig<-all_stats_OI[which(all_stats_OI$fdr<0.05 & abs(all_stats_OI$db)>=0.05),]
sig<-merge(sig, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")

betas_org<-as.data.frame(combat_organoid_Beta[which(rownames(combat_organoid_Beta)%in%CpG_OI),])
betas_org$CpG<-rownames(betas_org)
betas_org<-melt(betas_org)
betas_org$sample.type<-"organoid"
betas_org$sample.n<-"full"

betas_pri<-as.data.frame(ibd_combo_og[which(rownames(ibd_combo_og)%in%CpG_OI),])
betas_pri$CpG<-rownames(betas_pri)
betas_pri<-melt(betas_pri)
betas_pri$sample.type<-"primary"
betas_pri$sample.n<-"full"

betas_pri_epic<-as.data.frame(ibd_beta_epic[which(rownames(ibd_beta_epic)%in%CpG_OI),])
if(length(which(rownames(ibd_beta_epic)%in%CpG_OI))==1){
  colnames(betas_pri_epic)<-rownames(ibd_beta_epic)[which(rownames(ibd_beta_epic)%in%CpG_OI)]
  betas_pri_epic<-as.data.frame(t(betas_pri_epic))
  betas_pri_epic$CpG<-rownames(betas_pri_epic)
  betas_pri_epic<-melt(betas_pri_epic)
  betas_pri_epic$sample.type<-"primary"
  betas_pri_epic$sample.n<-"sub"
  
  betas<-rbind(betas_org,betas_pri,betas_pri_epic)}else{betas<-rbind(betas_org,betas_pri)}



betas<-merge(betas, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
plt<-merge(betas, meta_merge, by.x="variable", by.y="array.id")

plt<-plt[which(plt$sample.site=="TI"),]
###rm IBD-U from TI
plt<-plt[which(!(plt$variable%in%epic.organoid_exclude_other$array.id[which(epic.organoid_exclude_other$diagnosis=="IBD-U")])),]

sum_stat<-melt(tapply(plt$value, list(plt$MAPINFO,plt$diagnosis,plt$sample.site, plt$sample.type), mean))
colnames(sum_stat)<-c("MAPINFO","diagnosis","sample.site","sample.type","value")
sum_stat<-sum_stat[which(!(is.na(sum_stat$value))),]

colnames(sig)[6]<-"sample.site"


plt$new_color<-paste(plt$diagnosis, plt$sample.n)
plt<-plt[which(plt$CpG%in%sig$CpG),]

plt$diagnosis<-factor(plt$diagnosis, levels=c("Control","UC","CD"))
plt$CpG<-factor(plt$CpG, levels=CpG_OI)


p<-ggplot()+geom_boxplot(aes(diagnosis, value,fill=new_color, color=sample.n), plt, outlier.shape = NA)+facet_grid(sample.type~CpG)+
  geom_jitter(aes(diagnosis, value,fill=new_color, color=sample.n), plt, width=0.25, shape=21, size=2)+
  scale_color_manual(values=c("black","grey40"))+fillscale_diagnosis_sample+theme_bw()+th+ylim(0,1)+ylab("DNAm Beta Value")+xlab("Diagnosis")+
  theme(legend.position="none",plot.margin = margin(0.5, 0.1, 0.1, 0.5, "cm"),
        axis.text = element_text(size =12, color="black", hjust=1),
        axis.title = element_text(size =12),
        legend.text = element_text(size =12),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))

p
sig_plt<-plt[,c(2,4,7,8)][!duplicated(plt[,c(2,4,7,8)]),]
sig$sig<-"*"
sig_plt<-merge(sig_plt, sig[,c(1,5:7,9)],by=c("CpG","sample.type","diagnosis"),all.x=T)
sig_plt$CpG<-factor(sig_plt$CpG, levels=CpG_OI)

p+geom_text(aes(diagnosis,label = sig, y=0.95),data=sig_plt, size=8)

ggsave(file=here("DNAm/figs","celiac_TI_DNAm_box_threebatch.pdf"), width=5.5, height=7)
ggsave(file=here("DNAm/figs/jpeg","celiac_TI_DNAm_box_threebatch.jpeg"), width=7, height=5)



###########
## heatplot
###########
#' ## Load DNAm data
load(file=paste(here("DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))
load(here("DNAm/data/ibd_adjusted_combatted.Rdata"))
# load(file=paste(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))
# load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/ibd_adjusted_combatted.Rdata"))


## PCA together
identical(colnames(combat_ibd_Beta), sampleinfo_DNAm$array.id)
identical(colnames(combat_organoid_Beta), epic.organoid_combined$array.id)

combat_ibd_Beta<-combat_ibd_Beta[which(rownames(combat_ibd_Beta)%in%rownames(combat_organoid_Beta)),]
ibd_EPIC_organoid_beta_overlap<-combat_organoid_Beta[which(rownames(combat_organoid_Beta)%in%rownames(combat_ibd_Beta)),]
ibd_EPIC_organoid_beta_overlap<-ibd_EPIC_organoid_beta_overlap[match(rownames(combat_ibd_Beta),rownames(ibd_EPIC_organoid_beta_overlap)),]

identical(rownames(ibd_EPIC_organoid_beta_overlap),rownames(combat_ibd_Beta))
identical(colnames(combat_ibd_Beta), sampleinfo_DNAm$array.id)
identical(colnames(ibd_EPIC_organoid_beta_overlap), epic.organoid_combined$array.id)

ibd_both<-cbind(ibd_EPIC_organoid_beta_overlap,combat_ibd_Beta)
epic.organoid_combined$sample.type<-"organoid"
colnames(epic.organoid_combined)[which(colnames(epic.organoid_combined)=="passage.or.rescope.no_numeric")]<-"passage.or.rescope.no"
sampleinfo_both<-rbind(epic.organoid_combined[,c(2,1,3,4,18,12,5:7,9,10  )], sampleinfo_DNAm[,c(1,2,4:6,8, 11:13, 17,18)])
identical(colnames(ibd_both), sampleinfo_both$array.id)


pca<- prcomp(t(ibd_both[complete.cases(ibd_both),]))
vars <- pca$sdev^2
Importance<-vars/sum(vars)

Loadings<-as.data.frame(pca$x)
Loadings$array.id<-rownames(Loadings)

Loadings<-merge(Loadings, sampleinfo_both, by="array.id")

ggplot(Loadings, aes(PC1, PC2, fill=sample.site, color=sample.type))+geom_point(shape=21, size=3)+theme_bw()+
  fillscale_sampsite+
  th+scale_color_manual(values=c("white", "black")) + xlab("PC1 (25%)") + ylab("PC2 (22%)")

ggsave(here("DNAm/figs","PC1_PC2_Primary_organoid_threebatch.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_Primary_organoid_threebatch.jpeg"),  width = 7.5, height = 6)



#'### Primary
#' remove rescopes
sampleinfo_og<-sampleinfo_DNAm[which(is.na(sampleinfo_DNAm$passage.or.rescope.no)),]

MHCI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","MR1","CD1D","IRF1","NLRC5")
MHCI_sig_CpG<-unique(IBD_genes$CpG[which(IBD_genes$Gene.name%in%MHCI)])

MHCI_beta_sig<-ibd_both[which(rownames(ibd_both)%in%MHCI_sig_CpG),]

## exclued high passage IBD unknown and rescopes
print(table(epic.organoid_exclude_other$diagnosis))
epic.organoid_exclude_other$sample.type<-"organoid"
epic.organoid_exclude_other$passage.or.rescope.no<-as.character(epic.organoid_exclude_other$passage.or.rescope.no_numeric)


head(epic.organoid_exclude_other)
head(sampleinfo_DNAm[,c(1,2,4:6,8, 11:13, 17,18)])

sampleinfo_both<-rbind(epic.organoid_exclude_other[,c(1,14,15,2,28,29,16,17,3,5,6)], sampleinfo_DNAm[,c(1,2,4:6,8, 11:13, 17,18)])

keep<-c(sampleinfo_both$array.id[which(sampleinfo_both$sample.type=="organoid" & sampleinfo_both$array.id%in%epic.organoid_exclude_other$array.id)],
        sampleinfo_both$array.id[which(sampleinfo_both$sample.type=="purified" & sampleinfo_both$array.id%in%sampleinfo_og$array.id)])
sampleinfo_both<-sampleinfo_both[which(sampleinfo_both$array.id%in%keep),]

sampleinfo_both$diagnosis<-as.factor(sampleinfo_both$diagnosis)
print(levels(sampleinfo_both$diagnosis))
levels(sampleinfo_both$diagnosis)<-c("CD","Control","IBD-U","UC","CD","UC")

MHCI_beta_sig<-MHCI_beta_sig[,which(colnames(MHCI_beta_sig)%in%sampleinfo_both$array.id)]

heat_plot_df<-melt(MHCI_beta_sig)
heat_plot_df<-merge(heat_plot_df, sampleinfo_both, by.x="Var2", by.y="array.id")
heat_plot_df$Var2<-factor(heat_plot_df$Var2, levels=sampleinfo_both$array.id[order(sampleinfo_both$sample.type, sampleinfo_both$diagnosis)])

#' gene label
MHCI_genes<-IBD_genes[which(IBD_genes$Gene.name%in%MHCI),]
MHCI_genes<-MHCI_genes[,c("CpG","Gene.name")][!duplicated(MHCI_genes[,c("CpG","Gene.name")]),]
MHCI_genes<-as.data.frame(MHCI_genes%>%group_by(CpG) %>% summarize(label = paste(Gene.name, collapse = ',')))
MHCI_genes<-MHCI_genes[which(MHCI_genes$CpG%in%rownames(MHCI_beta_sig)),]
#MHCI_genes$CpG<-factor(MHCI_genes$CpG, levels=CpG_order)

CpG_order<-MHCI_genes$CpG[order(MHCI_genes$label)]

heat_plot_df$Var1<-factor(heat_plot_df$Var1, levels=CpG_order)
MHCI_genes$CpG<-factor(MHCI_genes$CpG, levels=CpG_order)

#beta_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))
beta_palette <- colorRampPalette((brewer.pal(9, "Blues")))


#heat_plot_df<-merge(heat_plot_df, MHCI_genes, by.x="Var1", by.y="CpG")


placeholder<-ggplot() + theme_void()

#margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")

myColors_type<- c("#bf812d","#542788")
color_possibilities_type<-c( "organoid","purified")
names(myColors_type) <- color_possibilities_type

fillscale_diagnosis_plus <- scale_fill_manual(name=NULL,values = c(myColors_diagnosis,myColors_type), drop = T)


heat_plt_MHC<-function(heat_plot_df, site){

  sampleinfo_both_TI<-sampleinfo_both[which(sampleinfo_both$sample.site==site),]
  sampleinfo_both_ordered<-sampleinfo_both_TI[order(sampleinfo_both_TI$sample.type, sampleinfo_both_TI$diagnosis),]
  org_CD_max<-max(which(sampleinfo_both_ordered$sample.type=="organoid" & sampleinfo_both_ordered$diagnosis=="CD"))
  pri_CD_max<-max(which(sampleinfo_both_ordered$sample.type=="purified" & sampleinfo_both_ordered$diagnosis=="CD"))
  pri_CD_min<-min(which(sampleinfo_both_ordered$sample.type=="purified" & sampleinfo_both_ordered$diagnosis=="CD"))

  squares<-ggplot()+geom_tile(aes(Var2, 1, fill=diagnosis),heat_plot_df)+
    geom_tile(aes(Var2, 2, fill=sample.type),heat_plot_df)+
    fillscale_diagnosis_plus+
    th+theme_void()+theme(plot.margin = unit(c(1, 0, 0, 0), "cm"),
                          legend.text = element_text(size = 7),
                          legend.position = "top",
                          legend.key.size = unit(0.3, 'cm'))


  sample_type_squares<-ggplot()+geom_tile(aes(Var2, 1, fill=sample.type),heat_plot_df)+
    th+theme_void()+scale_fill_manual(values=c("#bf812d","#542788"))+theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

  plot_grid(
    placeholder,squares,

    ggplot()+geom_text(aes(1, CpG, label=label),MHCI_genes, size=2.5,  hjust = 1)+xlim(0.7,1)+theme(plot.margin = unit(c(0, -1, 0, 0), "cm"))+
      theme(axis.text = element_blank(),axis.ticks = element_blank(), axis.title.y = element_blank(),
            panel.background = element_blank())+xlab(""),

    ggplot()+geom_tile(aes(Var2, Var1, fill=value),heat_plot_df)+
      ylab("CpG")+xlab("Individual")+theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(), axis.text.y = element_text(size=8))+
      scale_fill_gradientn(colours = beta_palette(100), limits=c(0, 1), name="DNAm\nBeta")+theme(plot.margin = unit(c(0, 0, 0, -1), "cm"))+
      geom_rect(aes(xmin=0,xmax=org_CD_max, ymin=0.5, ymax=length(CpG_order)+0.5), fill=NA, color="black")+
      geom_rect(aes(xmin=pri_CD_min,xmax=pri_CD_max, ymin=0.5, ymax=length(CpG_order)+0.5), fill=NA, color="black"),

    align="vh",axis="rl", nrow=2, rel_heights = c(1,6), rel_widths = c(1,2))
}

heat_plt_MHC(heat_plot_df[which(heat_plot_df$sample.site=="TI"),], "TI")
ggsave(file=here("DNAm/figs", "heat_plot_MHCI_TI_threebatch.pdf"),heat_plt_MHC(heat_plot_df[which(heat_plot_df$sample.site=="TI"),], "TI"), w=8, h=7 )
ggsave(file=here("DNAm/figs/jpeg", "heat_plot_MHCI_TI_threebatch.jpeg"),heat_plt_MHC(heat_plot_df[which(heat_plot_df$sample.site=="TI"),], "TI"), w=8, h=7 )

heat_plt_MHC(heat_plot_df[which(heat_plot_df$sample.site=="SC"),], "SC")
ggsave(file=here("DNAm/figs", "heat_plot_MHCI_SC_threebatch.pdf"),heat_plt_MHC(heat_plot_df[which(heat_plot_df$sample.site=="SC"),], "SC"), w=8, h=7 )
ggsave(file=here("DNAm/figs/jpeg", "heat_plot_MHCI_SC_threebatch.jpeg"),heat_plt_MHC(heat_plot_df[which(heat_plot_df$sample.site=="SC"),], "SC"), w=8, h=7 )

print("HEAT PLOTS UP TO DATE")


###################
#'### Manhattan
###################

#' ## delta beta correlation plot
manhattan_data<-function(CpG, pvalue, fdr, db, diagnosis, tissue, sample.type){
  data.frame(CpG=CpG, db=db, p.value=pvalue, fdr=fdr, diagnosis=diagnosis, tissue=tissue, sample.type=sample.type)
}

sig_hits<-function(manhattan_data_df){
  manhattan_data_df[which(manhattan_data_df$fdr<0.05 & abs(manhattan_data_df$db)>0.05),]
}

load(file=here("DNAm/data/diff_DNAm_TI_UCCTRL.RData"))
load(here("DNAm/data/diff_DNAm_TI_UCCTRL_organoid_ThreeBatch.RData"))

diff_meth_TI_UCCTRL$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL))
diff_meth_TI_UCCTRL_organoid$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL_organoid$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL_organoid))

manhattan_TI_CD_primary<-manhattan_data(rownames(diff_meth_TI_UCCTRL), diff_meth_TI_UCCTRL$p.value, diff_meth_TI_UCCTRL$adjusted_p, diff_meth_TI_UCCTRL$delta_beta, "CD", "TI", "primary")
manhattan_TI_CD_organoid<-manhattan_data(rownames(diff_meth_TI_UCCTRL_organoid), diff_meth_TI_UCCTRL_organoid$p.value, diff_meth_TI_UCCTRL_organoid$adjusted_p, diff_meth_TI_UCCTRL_organoid$delta_beta, "CD", "TI", "organoid")

## EPIC only CpGs
load(file=here("DNAm/data","diff_DNAm_EPIC_only.RData"))
diff_meth_TI_UCCTRL$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL))
manhattan_TI_CD_primary_EPIC<-manhattan_data(rownames(diff_meth_TI_UCCTRL), diff_meth_TI_UCCTRL$p.value, diff_meth_TI_UCCTRL$adjusted_p, diff_meth_TI_UCCTRL$delta_beta, "CD", "TI", "primary")

# downloaded from illumina website
anno_EPIC<-read.csv(here("DNAm/data","MethylationEPIC_v-1-0_B4.csv"), skip=7)


manhattan_TI_CD_primary<-merge(manhattan_TI_CD_primary, anno_EPIC[,c("IlmnID","Genome_Build","CHR","MAPINFO")], by.x="CpG", by.y="IlmnID")
manhattan_TI_CD_organoid<-merge(manhattan_TI_CD_organoid, anno_EPIC[,c("IlmnID","Genome_Build","CHR","MAPINFO")], by.x="CpG", by.y="IlmnID")
manhattan_TI_CD_primary_EPIC<-merge(manhattan_TI_CD_primary_EPIC, anno_EPIC[,c("IlmnID","Genome_Build","CHR","MAPINFO")], by.x="CpG", by.y="IlmnID")


      #manhattan_TI_CD_primary<-manhattan_TI_CD_primary[sample(1:nrow(manhattan_TI_CD_primary), 1000),]


## posistion cumsum
anno_EPIC_max_position<-as.data.frame(anno_EPIC %>% group_by(CHR) %>% slice(which.max(MAPINFO)))
anno_EPIC_max_position<-anno_EPIC_max_position[,c("IlmnID","Genome_Build","CHR","MAPINFO")]
anno_EPIC_max_position<-anno_EPIC_max_position[which(!(anno_EPIC_max_position$CHR%in%c("X","Y"))),]
anno_EPIC_max_position$CHR<-as.numeric(as.character(anno_EPIC_max_position$CHR))
anno_EPIC_max_position<-anno_EPIC_max_position[order(anno_EPIC_max_position$CHR),]
anno_EPIC_max_position$MAPINFO<-as.numeric(as.character(anno_EPIC_max_position$MAPINFO))
anno_EPIC_max_position$CumSumPosistion<-cumsum(anno_EPIC_max_position$MAPINFO)


manhattan_TI_CD_primary$cumulative_posistion<-sapply(1:nrow(manhattan_TI_CD_primary), function(x) {
  chr<-as.numeric(manhattan_TI_CD_primary$CHR[x])
  if(chr==1){add=0}else{add<-as.numeric(anno_EPIC_max_position[,c("CumSumPosistion")][which(anno_EPIC_max_position$CHR==(chr-1))])}
  manhattan_TI_CD_primary$MAPINFO[x]+add
})

manhattan_TI_CD_organoid$cumulative_posistion<-sapply(1:nrow(manhattan_TI_CD_organoid), function(x) {
  chr<-as.numeric(manhattan_TI_CD_organoid$CHR[x])
  if(chr==1){add=0}else{add<-as.numeric(anno_EPIC_max_position[,c("CumSumPosistion")][which(anno_EPIC_max_position$CHR==(chr-1))])}
  manhattan_TI_CD_organoid$MAPINFO[x]+add
})


## point color
TI_CD_primary_sig<-sig_hits(manhattan_TI_CD_primary)
TI_CD_organoid_sig<-sig_hits(manhattan_TI_CD_organoid)

manhattan_TI_CD_primary$color<-sapply(1:nrow(manhattan_TI_CD_primary), function(x) {
  if((as.numeric(manhattan_TI_CD_primary$CHR[x]) %% 2) == 0){"even"}else{"odd"}})

manhattan_TI_CD_primary$color<-sapply(1:nrow(manhattan_TI_CD_primary), function(x) {
  if(manhattan_TI_CD_primary$fdr[x]<=0.05 & abs(manhattan_TI_CD_primary$db[x])>0.05){
    if(manhattan_TI_CD_primary$db[x]>=0){"Less DNAm in CD"}else{"More DNAm in CD"}
  }else{manhattan_TI_CD_primary$color[x]}})

manhattan_TI_CD_organoid$color<-sapply(1:nrow(manhattan_TI_CD_organoid), function(x) {
  if((as.numeric(manhattan_TI_CD_organoid$CHR[x]) %% 2) == 0){"even"}else{"odd"}})

manhattan_TI_CD_organoid$color<-sapply(1:nrow(manhattan_TI_CD_organoid), function(x) {
  if(manhattan_TI_CD_organoid$fdr[x]<=0.05 & abs(manhattan_TI_CD_organoid$db[x])>0.05){
    if(manhattan_TI_CD_organoid$db[x]>=0){"Less DNAm in CD"}else{"More DNAm in CD"}
  }else{manhattan_TI_CD_organoid$color[x]}})

anno_EPIC_max_position$axis_pos<-(anno_EPIC_max_position$MAPINFO/2)+anno_EPIC_max_position$CumSumPosistion

anno_EPIC_max_position$axis_pos<-sapply(1:nrow(anno_EPIC_max_position), function(x){
  chr<-anno_EPIC_max_position$CHR[x]
  if(chr==1){(anno_EPIC_max_position$MAPINFO[x]/2)}else{
    (anno_EPIC_max_position$MAPINFO[x-1]/2)+anno_EPIC_max_position$CumSumPosistion[x-1]}
})

ggplot(manhattan_TI_CD_primary, aes(cumulative_posistion,-log10(p.value)))+geom_point(aes(color=color), size=0.75)+
  geom_hline(yintercept = -log10(1.984530e-05), color="grey")+
  theme_bw()+th_present+scale_color_manual(values=c("grey60","blue","red","grey80"), name="Differential\nDNAm\nDirection", guide=F)+
  xlab("Chromosome")+ylab("P Value (-log10)")+
  scale_x_continuous(label = anno_EPIC_max_position$CHR, breaks = anno_EPIC_max_position$axis_pos)

ggsave(file=here("DNAm/figs/jpeg", "Manhattan_primary_threebatch.jpeg"),w=15, h=4 )

ggplot(manhattan_TI_CD_organoid, aes(cumulative_posistion,-log10(p.value)))+geom_point(aes(color=color), size=0.75)+
  geom_hline(yintercept = -log10(0.0001667903), color="grey")+
  theme_bw()+th_present+scale_color_manual(values=c("grey60","blue","red","grey80"), name="Differential\nDNAm\nDirection", guide=F)+
  xlab("Chromosome")+ylab("P Value (-log10)")+
  scale_x_continuous(label = anno_EPIC_max_position$CHR, breaks = anno_EPIC_max_position$axis_pos)

ggsave(file=here("DNAm/figs/jpeg", "Manhattan_organoid_threebatch.jpeg"),w=15, h=4 )

plot_grid(ggplot(manhattan_TI_CD_organoid, aes(cumulative_posistion,-log10(p.value)))+geom_point(aes(color=color), size=0.75)+
            geom_hline(yintercept = -log10(0.0001667903), color="grey")+
            theme_bw()+th_present+scale_color_manual(values=c("grey60","blue","red","grey80"), name="Differential\nDNAm\nDirection", guide=F)+
            xlab("Chromosome")+ylab("P Value (-log10)")+
            scale_x_continuous(label = anno_EPIC_max_position$CHR, breaks = anno_EPIC_max_position$axis_pos),
          ggplot(manhattan_TI_CD_primary, aes(cumulative_posistion,-log10(p.value)))+geom_point(aes(color=color), size=0.75)+
            geom_hline(yintercept = -log10(1.984530e-05), color="grey")+
            theme_bw()+th_present+scale_color_manual(values=c("grey60","blue","red","grey80"), name="Differential\nDNAm\nDirection", guide=F)+
            xlab("Chromosome")+ylab("P Value (-log10)")+
            scale_x_continuous(label = anno_EPIC_max_position$CHR, breaks = anno_EPIC_max_position$axis_pos),
          align="h", ncol=1)

ggsave(file=here("DNAm/figs/jpeg", "Manhattan_both_aligned_threebatch.jpeg"),w=15, h=8 )

#############
### Manhattan MHC region
#############
MHCI_region_primary<-manhattan_TI_CD_primary[which(manhattan_TI_CD_primary$CHR=="6"),]
MHCI_region_primary$MAPINFO<-as.numeric(MHCI_region_primary$MAPINFO)
MHCI_region_primary<-MHCI_region_primary[which(MHCI_region_primary$MAPINFO>=29670000 & MHCI_region_primary$MAPINFO<=33000000),]
dim(MHCI_region_primary)
MHCI_region_organoid<-manhattan_TI_CD_organoid[which(manhattan_TI_CD_organoid$CHR=="6"),]
MHCI_region_organoid$MAPINFO<-as.numeric(MHCI_region_organoid$MAPINFO)
MHCI_region_organoid<-MHCI_region_organoid[which(MHCI_region_organoid$MAPINFO>=29670000 & MHCI_region_organoid$MAPINFO<=33000000),]
dim(MHCI_region_organoid)



plot_grid(ggplot(MHCI_region_organoid, aes(MAPINFO/1000000,-log10(p.value)))+geom_point(aes(color=color), size=0.75)+
            geom_hline(yintercept = -log10(0.0001667903), color="grey")+
            theme_bw()+th_present+scale_color_manual(values=c("grey60","blue","red","grey80"), name="Differential\nDNAm\nDirection", guide=F)+
            xlab("Genomic Coordinate (Mb)")+ylab("P Value (-log10)"),
          ggplot(MHCI_region_primary, aes(MAPINFO/1000000,-log10(p.value)))+geom_point(aes(color=color), size=0.75)+
            geom_hline(yintercept = -log10(1.984530e-05), color="grey")+
            theme_bw()+th_present+scale_color_manual(values=c("grey60","blue","red","grey80"), name="Differential\nDNAm\nDirection", guide=F)+
            xlab("Genomic Coordinate (Mb)")+ylab("P Value (-log10)"),
          align="h", ncol=1)

ggsave(file=here("DNAm/figs/jpeg", "Manhattan_both_aligned_zoom_threebatch.jpeg"),w=10, h=8 )


###################
#'### Where are the CpGs
###################
EPIC_featureCount<-tapply(EPIC_genes$CpG, EPIC_genes$CpG_in, length)
EPIC_featureCount<-data.frame(CpG_Count=as.numeric(EPIC_featureCount),
                              Feature=names(EPIC_featureCount))
EPIC_featureCount$Feature<-factor(EPIC_featureCount$Feature,
                                  levels=c("promoter","intragenic","3'UTR"))
levels(EPIC_featureCount$Feature)<-c("Promoter","Gene Body","3'UTR")

EPIC_featureCount$percent<-(EPIC_featureCount$CpG_Count/sum(EPIC_featureCount$CpG_Count))*100

ggplot(EPIC_featureCount, aes(Feature, CpG_Count, fill=Feature))+
  geom_bar(position=position_dodge(width=0.9),stat="identity", color="grey25")+theme_bw()+
  ylab("CpG Count")+xlab("")+th+
  theme(legend.position="none",
        axis.text = element_text(size =10, color="black"),
        axis.title = element_text(size =15),
        legend.text = element_text(size =14),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))+
  geom_text(aes(label = paste (round(percent, 0), "%"), y=CpG_Count, vjust=-0.5), size=4, color="grey30")+
  scale_fill_manual(values=c("#d8d8d8","#cb95cc","#44dfcf","#fac900","#ff6969","#ff0000","grey97"))+
  scale_y_continuous(label=comma)


#'## TI MHC
sig_CpG_TICD<-rbind(sig_hits(manhattan_TI_CD_primary),sig_hits(manhattan_TI_CD_organoid),
                    sig_hits(manhattan_TI_CD_primary_EPIC))

IBD_genes_TICD<-merge(sig_CpG_TICD, EPIC_genes[,c(3,6,7,9)], by.x="CpG", by.y="IlmnID")
IBD_genes_TICD_MHC<-IBD_genes_TICD[which(IBD_genes_TICD$Gene.name%in%MHCI),c("CpG","Gene.name","CpG_in")]
IBD_genes_TICD_MHC<-IBD_genes_TICD_MHC[!(duplicated(IBD_genes_TICD_MHC)),]

length(unique(IBD_genes_TICD_MHC$CpG))


MHC_featureCount<-tapply(IBD_genes_TICD_MHC$CpG, IBD_genes_TICD_MHC$CpG_in, length)
MHC_featureCount<-data.frame(CpG_Count=as.numeric(MHC_featureCount),
                              Feature=names(MHC_featureCount))
MHC_featureCount$Feature<-factor(MHC_featureCount$Feature,
                                  levels=c("promoter","intragenic","3'UTR"))
levels(MHC_featureCount$Feature)<-c("Promoter","Gene Body","3'UTR")

MHC_featureCount$percent<-(MHC_featureCount$CpG_Count/sum(MHC_featureCount$CpG_Count))*100

ggplot(MHC_featureCount, aes(Feature, CpG_Count, fill=Feature))+
  geom_bar(position=position_dodge(width=0.9),stat="identity", color="grey25")+theme_bw()+
  ylab("CpG Count")+xlab("")+th+
  theme(legend.position="none",
        axis.text = element_text(size =10, color="black"),
        axis.title = element_text(size =15),
        legend.text = element_text(size =14),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))+
  geom_text(aes(label = paste (round(percent, 0), "%"), y=CpG_Count, vjust=-0.5), size=4, color="grey30")+
  scale_fill_manual(values=c("#d8d8d8","#cb95cc","#44dfcf","#fac900","#ff6969","#ff0000","grey97"))+
  scale_y_continuous(label=comma)

#'## R Session Info
sessionInfo()


