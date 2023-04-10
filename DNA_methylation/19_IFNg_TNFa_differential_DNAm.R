#'---
#'title: Differentiation and treatment organoids
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
#suppressMessages(library(minfi))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(RColorBrewer))
suppressMessages(library(here))
suppressMessages(library(binom))
suppressMessages(library(limma))
# suppressMessages(library(sva))
# suppressMessages(library(pamr))
# suppressMessages(library(GEOquery))
# suppressMessages(library(GEOmetadb))

suppressMessages(library(dplyr))
#suppressMessages(library(lmtest))
suppressMessages(library(gridExtra))
suppressMessages(library(gtools))
#suppressMessages(library(rafalib))
library(ggsignif)
library(cowplot)



options(stringsAsFactors = FALSE)


#' ### Load Functions
source(here("general_functions","00_pretty_plots.R"))



#load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/desktop_june2022/Documents/DNAm_organoid_passage/data/validation/DNAm","validation_betas_normalized.RData"))
load(file=here("../../DNAm_organoid_passage/data/validation/DNAm","validation_betas_normalized.RData"))


#' 
#' ################################################################################
#' #' ## Differential methylation with differentiation 
#' ################################################################################
#' #########
#' #' ## low passage
#' #########
#' validation_epic.organoid_low<-validation_epic.organoid[which(validation_epic.organoid$passage_hilo=="low" & validation_epic.organoid$comparison=="differentiation"),]
#' table(validation_epic.organoid_low$differentiation, validation_epic.organoid_low$individual)
#' 
#' validation_organoid_beta_low<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_low$array.id)]
#' identical(validation_epic.organoid_low$array.id, colnames(validation_organoid_beta_low))
#' 
#' validation_epic.organoid_low$individual<-as.factor(validation_epic.organoid_low$individual)
#' 
#' mod<-model.matrix(~ differentiation + individual, data=validation_epic.organoid_low)
#' fit <- lmFit(validation_organoid_beta_low, mod)
#' ebfit <- eBayes(fit)
#' 
#' # Delta beta
#' low_pass_diff_db<-sapply(1:nrow(validation_organoid_beta_low), function(x){
#'   sampleinfo_cpg<-validation_epic.organoid_low
#'   sampleinfo_cpg$beta<-as.numeric(validation_organoid_beta_low[x,])
#'   mns<-tapply(sampleinfo_cpg$beta,  sampleinfo_cpg$differentiation, mean)
#'   mns[2]-mns[1]# UD minus D
#' })
#' 
#' 
#' differnetiation_low_stats<-data.frame(p.value=ebfit$p.value[,"differentiationUD"], CpG=rownames(validation_organoid_beta_low), db=low_pass_diff_db)
#' 
#' # Adjust P values
#' differnetiation_low_stats$p_adjusted<-p.adjust(differnetiation_low_stats$p.value, method="BH")
#' 
#' diff_CpG_lowdiff<-differnetiation_low_stats[which(differnetiation_low_stats$p_adjusted<0.05 & abs(differnetiation_low_stats$db)>0.15),] #
#' diff_CpG_db_lowdiff<-diff_CpG_lowdiff$CpG[which((diff_CpG_lowdiff$db)>=0.15)] #  
#' diff_CpG_db_lowdiff<-diff_CpG_lowdiff$CpG[which((diff_CpG_lowdiff$db)<=(-0.15))] #  
#' 
#' 
#' 
#' #########
#' #' ## high passage
#' #########
#' validation_epic.organoid_high<-validation_epic.organoid[which(validation_epic.organoid$passage_hilo=="high" & validation_epic.organoid$comparison=="differentiation"),]
#' table(validation_epic.organoid_high$differentiation, validation_epic.organoid_high$individual)
#' 
#' validation_organoid_beta_high<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_high$array.id)]
#' identical(validation_epic.organoid_high$array.id, colnames(validation_organoid_beta_high))
#' 
#' validation_epic.organoid_high$individual<-as.factor(validation_epic.organoid_high$individual)
#' 
#' mod<-model.matrix(~ differentiation  + individual, data=validation_epic.organoid_high)
#' fit <- lmFit(validation_organoid_beta_high, mod)
#' ebfit <- eBayes(fit)
#' 
#' # Delta beta
#' high_pass_diff_db<-sapply(1:nrow(validation_organoid_beta_high), function(x){
#'   sampleinfo_cpg<-validation_epic.organoid_high
#'   sampleinfo_cpg$beta<-as.numeric(validation_organoid_beta_high[x,])
#'   mns<-tapply(sampleinfo_cpg$beta,  sampleinfo_cpg$differentiation, mean)
#'   mns[2]-mns[1]# UD minus D
#' })
#' 
#' 
#' differnetiation_high_stats<-data.frame(p.value=ebfit$p.value[,"differentiationUD"], CpG=rownames(validation_organoid_beta_high), db=high_pass_diff_db)
#' 
#' # Adjust P values
#' differnetiation_high_stats$p_adjusted<-p.adjust(differnetiation_high_stats$p.value, method="BH")
#' 
#' diff_CpG_highdiff<-differnetiation_high_stats[which(differnetiation_high_stats$p_adjusted<0.05 & abs(differnetiation_high_stats$db)>0.15),] #
#' diff_CpG_db_highdiff<-diff_CpG_highdiff$CpG[which((diff_CpG_highdiff$db)>=0.15)] #  
#' diff_CpG_db_highdiff<-diff_CpG_highdiff$CpG[which((diff_CpG_highdiff$db)<=(-0.15))] #  
#' 
#' 
#' 
#' #' 
#' #' #' ### CpG to gene associations
#' #' EPIC_genes<-read.csv(here("data","EPIC_ensembl_gene_annotation.csv")) # 1137194
#' #' 
#' #' diff_genes_db_hypovalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypovalidation)] ) #11442
#' #' diff_genes_db_hypervalidation<-unique(EPIC_genes$Gene.name[which(EPIC_genes$IlmnID%in%diff_CpG_db_hypervalidation)] ) # 2084
#' #' 
#' #' write.table(diff_genes_db_hypovalidation, file=here("data/validation/DNAm/","validation_genes_hypomethylation.txt"), quote=F, row.names = F, col.names = F)
#' #' write.table(diff_genes_db_hypervalidation, file=here("data/validation/DNAm/","validation_genes_hypermethylation.txt"), quote=F, row.names = F, col.names = F)
#' 


#######################################################################
#' ## Differential methylation with IFNg 
#######################################################################
#' ## low passage
validation_epic.organoid_IFNg<-validation_epic.organoid[which(validation_epic.organoid$comparison=="cytokine" & validation_epic.organoid$treatment!="TNFa"),]
table(validation_epic.organoid_IFNg$treatment, validation_epic.organoid_IFNg$individual)

validation_organoid_beta_IFNg<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_IFNg$array.id)]
identical(validation_epic.organoid_IFNg$array.id, colnames(validation_organoid_beta_IFNg))

validation_epic.organoid_IFNg$individual<-as.factor(validation_epic.organoid_IFNg$individual)

mod<-model.matrix(~ treatment  + individual, data=validation_epic.organoid_IFNg)
fit <- lmFit(validation_organoid_beta_IFNg, mod)
ebfit <- eBayes(fit)

# Delta beta
low_pass_IFNg_db<-sapply(1:nrow(validation_organoid_beta_IFNg), function(x){
  sampleinfo_cpg<-validation_epic.organoid_IFNg
  sampleinfo_cpg$beta<-as.numeric(validation_organoid_beta_IFNg[x,])
  mns<-tapply(sampleinfo_cpg$beta,  sampleinfo_cpg$treatment, mean)
  mns[2]-mns[1]# UT minus IFNg
})


IFNg_stats<-data.frame(p.value=ebfit$p.value[,"treatmentUT"], CpG=rownames(validation_organoid_beta_IFNg), db=low_pass_IFNg_db)

# Adjust P values
IFNg_stats$p_adjusted<-p.adjust(IFNg_stats$p.value, method="BH")

diff_CpG_IFNg<-IFNg_stats[which(IFNg_stats$p_adjusted<0.05 & abs(IFNg_stats$db)>0.15),] #
diff_CpG_IFNg_hypo<-diff_CpG_IFNg$CpG[which((diff_CpG_IFNg$db)>=0.15)] #  
diff_CpG_IFNg_hyper<-diff_CpG_IFNg$CpG[which((diff_CpG_IFNg$db)<=(-0.15))] #  


#######################################################################
#' ## Differential methylation with TNFa 
#######################################################################
#' ## low passage
validation_epic.organoid_TNFa<-validation_epic.organoid[which( validation_epic.organoid$comparison=="cytokine" & validation_epic.organoid$treatment!="IFNg"),]
table(validation_epic.organoid_TNFa$treatment, validation_epic.organoid_TNFa$individual)

validation_organoid_beta_TNFa<-validation_organoid_beta[,which(colnames(validation_organoid_beta)%in%validation_epic.organoid_TNFa$array.id)]
identical(validation_epic.organoid_TNFa$array.id, colnames(validation_organoid_beta_TNFa))

validation_epic.organoid_TNFa$individual<-as.factor(validation_epic.organoid_TNFa$individual)

mod<-model.matrix(~ treatment  + individual, data=validation_epic.organoid_TNFa)
fit <- lmFit(validation_organoid_beta_TNFa, mod)
ebfit <- eBayes(fit)

# Delta beta
low_pass_TNFa_db<-sapply(1:nrow(validation_organoid_beta_TNFa), function(x){
  sampleinfo_cpg<-validation_epic.organoid_TNFa
  sampleinfo_cpg$beta<-as.numeric(validation_organoid_beta_TNFa[x,])
  mns<-tapply(sampleinfo_cpg$beta,  sampleinfo_cpg$treatment, mean)
  mns[2]-mns[1]# UT minus TNFa
})


TNFa_stats<-data.frame(p.value=ebfit$p.value[,"treatmentUT"], CpG=rownames(validation_organoid_beta_TNFa), db=low_pass_TNFa_db)

# Adjust P values
TNFa_stats$p_adjusted<-p.adjust(TNFa_stats$p.value, method="BH")

diff_CpG_TNFa<-TNFa_stats[which(TNFa_stats$p_adjusted<0.05 & abs(TNFa_stats$db)>0.15),] #
diff_CpG_TNFa_hypo<-diff_CpG_TNFa$CpG[which((diff_CpG_TNFa$db)>=0.15)] #  
diff_CpG_TNFa_hyper<-diff_CpG_TNFa$CpG[which((diff_CpG_TNFa$db)<=(-0.15))] #  



### plot whatever gene DMR
EPIC_genes<-read.csv(here("DNAm/data","EPIC_ensembl_gene_annotation.csv")) # 1137194
#EPIC_genes<-read.csv(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","EPIC_ensembl_gene_annotation.csv")) # 1137194

gene_DNAm<-function(gene){
  CpG_goi<-EPIC_genes[which(EPIC_genes$Gene.name%in%gene & EPIC_genes$CpG_in=="promoter"),]
  CpGs<-unique(CpG_goi$IlmnID)

  meta_TI<-validation_epic.organoid[which(validation_epic.organoid$comparison=="cytokine" ),]

  betas<-validation_organoid_beta[which(rownames(validation_organoid_beta)%in%CpGs),]
  betas<-melt(betas)
  organoid_plt<-merge(meta_TI, betas, by.x="array.id",by.y="X2")
  organoid_plt<-merge(organoid_plt, CpG_goi, by.x="X1",by.y="IlmnID")
  
  organoid_plt$treatment<-factor(organoid_plt$treatment, levels=c("UT","IFNg","TNFa"))
  
  ggplot(organoid_plt, aes(X1, value, fill=treatment))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(size=1,shape=21, color="black",position = position_dodge(width=0.75))+
    th+theme_bw()+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
    scale_fill_manual(values=c("grey80","cornflowerblue","firebrick4"), name="Treatment")+
    ylab("DNAm Beta")+xlab("")+ylim(0,1)}


gene_DNAm("NLRC5")

ggsave(here("DNAm/figs","DNAm_NLRC5_treatments_organoids.pdf"),width = 10, height = 5)
ggsave(here("DNAm/figs/jpeg","DNAm_NLRC5_treatments_organoids.jpeg"), width = 10, height = 5)








CpGs<-c("cg07839457","cg07862320")
meta_TI<-validation_epic.organoid[which(validation_epic.organoid$comparison=="cytokine" ),]
betas<-validation_organoid_beta[which(rownames(validation_organoid_beta)%in%CpGs),]
betas<-melt(betas)
organoid_plt<-merge(meta_TI, betas, by.x="array.id",by.y="X2")
organoid_plt$treatment<-factor(organoid_plt$treatment, levels=c("UT","IFNg","TNFa"))

comp<-list(c("IFNg","UT"),c("TNFa","IFNg"),c("TNFa","UT"))
ggplot(organoid_plt, aes(treatment, value, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size=2,shape=21, color="black",position = position_dodge(width=0.75))+
  th+theme_bw()+
  scale_fill_manual(values=c("grey80","cornflowerblue","firebrick4"), name="Treatment")+
  ylab("DNAm Beta")+xlab("")+ylim(0,1)+facet_wrap(~X1)+
  theme(legend.position="none",plot.margin = margin(0.5, 0.1, 0.1, 0.5, "cm"),
        axis.text = element_text(size =10, angle=45,color="black", hjust=1),
        axis.title = element_text(size =12),
        legend.text = element_text(size =12),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))+
  geom_signif(comparisons = comp, step_increase = 0.05,tip_length = 0.01,
              size = 0.5,vjust = 0.6,
              textsize = 3,  map_signif_level = T, color="grey60")




ggsave(here("DNAm/figs","DNAm_NLRC5_treatments_organoids_mainCpGs.pdf"),width = 5, height = 2.5)
ggsave(here("DNAm/figs/jpeg","DNAm_NLRC5_treatments_organoids_mainCpGs.jpeg"), width = 5, height = 2.5)




#####################
## plot all CpGs shown in paper
#####################


CpG_treatment_plot<-function(CpGs, gene){
  meta_TI<-validation_epic.organoid[which(validation_epic.organoid$comparison=="cytokine" ),]
  betas<-validation_organoid_beta[which(rownames(validation_organoid_beta)%in%CpGs),]
  betas<-melt(betas)
  organoid_plt<-merge(meta_TI, betas, by.x="array.id",by.y="X2")
  organoid_plt$treatment<-factor(organoid_plt$treatment, levels=c("UT","IFNg","TNFa"))
  
  comp<-list(c("IFNg","UT"),c("TNFa","IFNg"),c("TNFa","UT"))
  ggplot(organoid_plt, aes(treatment, value, fill=treatment))+
    geom_boxplot(outlier.shape = NA)+
    geom_point(size=2,shape=21, color="black",position = position_dodge(width=0.75))+
    th+theme_bw()+
    scale_fill_manual(values=c("grey80","cornflowerblue","firebrick4"), name="Treatment")+
    ylab("DNA Methylation")+xlab("")+ylim(0,1)+facet_wrap(~X1)+
    theme(legend.position="none",
          axis.text = element_text(size =10),
          axis.title = element_text(size =12),
          legend.text = element_text(size =12),
          legend.title = element_text(size =12),
          strip.text.x = element_text(size = 12))+
    geom_signif(comparisons = comp, step_increase = 0.05,tip_length = 0.01,
                size = 0.5,vjust = 0.6,
                textsize = 3,  map_signif_level = T, color="grey60")+
    ggtitle(gene)}


save_p<-plot_grid(CpG_treatment_plot(c("cg07839457","cg07862320"), "NLRC5"),
          CpG_treatment_plot(c("cg11706729","cg02756056"), "TAP1"),
          CpG_treatment_plot(c("cg23235965","cg11594821"), "HLA-E"),
          CpG_treatment_plot(c("cg27537252","cg05475649"), "B2M"))

ggsave(file=here("DNAm/figs/","keyMHCI_Treated.pdf"),save_p, width=6, height=6)



save_p<-plot_grid(
  CpG_treatment_plot(c("cg23533285","cg06422189"), "HLA-B"),
  CpG_treatment_plot(c("cg01064373","cg18511546"), "HLA-C"),
  CpG_treatment_plot(c("cg01405582","cg11587584"), "HLA-F"),
  CpG_treatment_plot(c("cg15375424","cg11666365"), "IRF1"),
  CpG_treatment_plot(c("cg24898914","cg16890093"), "PSMB8"),
  CpG_treatment_plot(c("cg06791592","cg10817441"), "PSMB9"),
  CpG_treatment_plot(c("cg18243910","cg27618777"), "TAP2"),
  CpG_treatment_plot(c("cg22834109","cg04392234"), "IL8"),ncol=2)
save_p

ggsave(file=here("DNAm/figs/","keyMHCI_supplement_Treated.pdf"),save_p, width=6, height=13)

CpG_treatment_plot(EPIC_genes[grep("IL8", EPIC_genes$Gene.name),]$IlmnID, "IL8")



###################
### Scores from Fliss
###################
load("rna_seq/data/RE_bulkRNASeq_and_EPIC_Datasets.RData")


load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/desktop_june2022/Documents/DNAm_organoid_passage/data/validation/DNAm","validation_betas_normalized.RData"))
meta_TI<-validation_epic.organoid[which(validation_epic.organoid$comparison=="cytokine" ),]

re.dnam<-re.dnam[which(re.dnam$array.id%in%meta_TI$array.id),]

re.dnam$treatment.1<-as.factor(re.dnam$treatment.1)
levels(re.dnam$treatment.1)<-c("IFNg", "UT" ,  "TNFa")
re.dnam$treatment.1<-factor(re.dnam$treatment.1,levels=c( "UT" , "IFNg", "TNFa"))
                            
comp<-list(c("IFNg","UT"),c("TNFa","IFNg"),c("TNFa","UT"))
ggplot(re.dnam, aes(treatment.1, Avg.MHCI, fill=treatment.1))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size=2,shape=21, color="black",position = position_dodge(width=0.75))+
  th+theme_bw()+
  scale_fill_manual(values=c("grey80","cornflowerblue","firebrick4"), name="Treatment")+
  ylab("Average MHC I\nDNA Methylation")+xlab("")+
  theme(legend.position="none",
        axis.text = element_text(size =10),
        axis.title = element_text(size =12),
        legend.text = element_text(size =12),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))+
  geom_signif(comparisons = comp, step_increase = 0.05,tip_length = 0.01,
              size = 0.5,vjust = 0.6,
              textsize = 3,  map_signif_level = T, color="grey60")
ggsave(file=here("DNAm/figs/","averageMHCIDNAm_Treated.pdf"),width=3, height=4)





