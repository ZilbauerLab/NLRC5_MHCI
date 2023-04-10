#' ### Load Libraries
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(RColorBrewer))
suppressMessages(library(here))

suppressMessages(library(dplyr))
suppressMessages(library(gridExtra))
library(ggsignif)
library(cowplot)
library(DESeq2)


options(stringsAsFactors = FALSE)

#' ### Load Functions
source(here("general_functions","00_pretty_plots.R"))


###################
### counts from Fliss
###################
load("rna_seq/data/RE_GRCh38_RawCounts.RData")

### expression matched to n=18 samples in design and DNAm
pheno.re.ti_treatment<-pheno.re.ti[-grep("UD",pheno.re.ti$rnaseq.id),]


### IFNg
pheno.re.ti_IFNg<-pheno.re.ti_treatment[which(pheno.re.ti_treatment$treatment!="TNFa"),]
table(pheno.re.ti_IFNg$treatment)

raw.re.ti_IFNg<-raw.re.ti[,which(colnames(raw.re.ti)%in%pheno.re.ti_IFNg$rnaseq.id)]
identical(colnames(raw.re.ti_IFNg),pheno.re.ti_IFNg$rnaseq.id)


dds <- DESeqDataSetFromMatrix(countData = raw.re.ti_IFNg,
                              colData = pheno.re.ti_IFNg,
                              design = ~treatment)
dds

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$pvalue),]
sum(res$padj < 0.1, na.rm=TRUE)

resSig_IFNg <- res[which(res$padj<0.001 & abs(res$log2FoldChange) >1),]



### TNFa
pheno.re.ti_TNFa<-pheno.re.ti_treatment[which(pheno.re.ti_treatment$treatment!="IFNg"),]
table(pheno.re.ti_TNFa$treatment)

raw.re.ti_TNFa<-raw.re.ti[,which(colnames(raw.re.ti)%in%pheno.re.ti_TNFa$rnaseq.id)]
identical(colnames(raw.re.ti_TNFa),pheno.re.ti_TNFa$rnaseq.id)


dds <- DESeqDataSetFromMatrix(countData = raw.re.ti_TNFa,
                              colData = pheno.re.ti_TNFa,
                              design = ~treatment)
dds

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$pvalue),]
sum(res$padj < 0.1, na.rm=TRUE)

resSig_TNFa <- res[which(res$padj<0.001 & abs(res$log2FoldChange) >1),]


### TNFa vs IFNg
pheno.re.ti_TNFaIFNg<-pheno.re.ti_treatment[which(pheno.re.ti_treatment$treatment!="NT"),]
table(pheno.re.ti_TNFaIFNg$treatment)

raw.re.ti_TNFaIFNg<-raw.re.ti[,which(colnames(raw.re.ti)%in%pheno.re.ti_TNFaIFNg$rnaseq.id)]
identical(colnames(raw.re.ti_TNFaIFNg),pheno.re.ti_TNFaIFNg$rnaseq.id)


dds <- DESeqDataSetFromMatrix(countData = raw.re.ti_TNFaIFNg,
                              colData = pheno.re.ti_TNFaIFNg,
                              design = ~treatment)
dds

dds <- DESeq(dds)
res <- results(dds)
res

resOrdered <- res[order(res$pvalue),]
sum(res$padj < 0.1, na.rm=TRUE)

resSig_TNFaIFNg <- res[which(res$padj<0.001 & abs(res$log2FoldChange) >1),]



###########
#'# aggregate to gene level
###########
mart <- biomaRt::useMart(biomart = "ensembl",
                         dataset = "hsapiens_gene_ensembl")

ttg <- biomaRt::getBM(
  attributes = c("ensembl_transcript_id", "transcript_version",
                 "ensembl_gene_id", "external_gene_name", "description",
                 "transcript_biotype"),
  mart = mart)
ttg <- dplyr::rename(ttg, target_id = ensembl_transcript_id,
                     ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
ttg <- dplyr::select(ttg, c('target_id', 'ens_gene', 'ext_gene'))
head(ttg)

ttg_simple<-ttg[,c("ens_gene","ext_gene")]
ttg_simple<-ttg_simple[!duplicated(ttg_simple),]



###########
## sig genes
###########
resSig_IFNg<-as.data.frame(resSig_IFNg)
resSig_IFNg$ens_gene<-rownames(resSig_IFNg)
resSig_IFNg<-merge(resSig_IFNg, ttg_simple, by="ens_gene")

resSig_TNFa<-as.data.frame(resSig_TNFa)
resSig_TNFa$ens_gene<-rownames(resSig_TNFa)
resSig_TNFa<-merge(resSig_TNFa, ttg_simple, by="ens_gene")

resSig_TNFaIFNg<-as.data.frame(resSig_TNFaIFNg)
resSig_TNFaIFNg$ens_gene<-rownames(resSig_TNFaIFNg)
resSig_TNFaIFNg<-merge(resSig_TNFaIFNg, ttg_simple, by="ens_gene")


GOI = c('HLA-F', 'HLA-G', 'HLA-A', 'HLA-E', 'HLA-C', 'HLA-B',"TAP1","TAP2","PSMB9","PSMB8","B2M","IRF1","NLRC5","CXCL8")
resSig_IFNg[which(resSig_IFNg$ext_gene%in%GOI),]
resSig_TNFa[which(resSig_TNFa$ext_gene%in%GOI),]
resSig_TNFaIFNg[which(resSig_TNFaIFNg$ext_gene%in%GOI),]



###########
## plot genes
###########
gene_exp_plot_treatment<-function(gene){
  goi<-melt(as.data.frame(raw.re.ti[which(rownames(raw.re.ti)%in%(unique(ttg$ens_gene[which(ttg$ext_gene==gene)]))),]))
  colnames(goi)<-c("rnaseq.id","Count")
  
  plt<-merge(pheno.re.ti_treatment,goi, by="rnaseq.id")
  
  plt$treatment<-as.factor(plt$treatment)
  levels(plt$treatment)<-c("IFNg", "UT",   "TNFa")
  plt$treatment<-factor(plt$treatment,levels=c( "UT" , "IFNg", "TNFa"))
  
  
  comp<-list(c("IFNg","UT"),c("TNFa","IFNg"),c("TNFa","UT"))
  ggplot(plt, aes(treatment, Count,fill=treatment))+geom_boxplot()+
    geom_point(aes(fill=treatment), shape=21, color="black", size=2)+theme_bw()+th+
    scale_fill_manual(values=c("grey80","cornflowerblue","firebrick4"), name="Treatment")+
    geom_signif(comparisons = comp, step_increase = 0.05,tip_length = 0.01,
                size = 0.5,vjust = 0.6,
                textsize = 3,  map_signif_level = c("*"=0.01), color="grey60")+
    theme(legend.position="none",
          axis.text = element_text(size =10),
          axis.title = element_text(size =12),
          legend.text = element_text(size =12),
          legend.title = element_text(size =12),
          strip.text.x = element_text(size = 12))+
    ggtitle(gene)
}


save_p<-plot_grid(plotlist = lapply(GOI, gene_exp_plot_treatment), ncol=3, align="vh")
ggsave(file=here("rna_seq/figs/","keyMHCI_counts_Treated.pdf"),save_p, width=6, height=15)




###################
### Scores from Fliss
###################
load("rna_seq/data/RE_bulkRNASeq_and_EPIC_Datasets.RData")

load(file=here("/media/redgar/Seagate Portable Drive/EBI_backup/desktop_june2022/Documents/DNAm_organoid_passage/data/validation/DNAm","validation_betas_normalized.RData"))
meta_TI<-validation_epic.organoid[which(validation_epic.organoid$comparison=="cytokine" ),]

### expression matched to n=18 samples in design and DNAm
pheno.re<-pheno.re[which(pheno.re$tissue=="TI"),]
pheno.re<-pheno.re[which(pheno.re$differentiated=="undifferentiated"),]
pheno.re<-pheno.re[which(pheno.re$sample.id%in%re.dnam$sample.id),]
pheno.re<-pheno.re[-grep("UD",pheno.re$rnaseq.id),]

pheno.re$treatment<-as.factor(pheno.re$treatment)
levels(pheno.re$treatment)<-c("IFNg", "UT" ,  "TNFa")
pheno.re$treatment<-factor(pheno.re$treatment,levels=c( "UT" , "IFNg", "TNFa"))

comp<-list(c("IFNg","UT"),c("TNFa","IFNg"),c("TNFa","UT"))
ggplot(pheno.re, aes(treatment, Avg.MHCI.exprs, fill=treatment))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size=2,shape=21, color="black",position = position_dodge(width=0.75))+
  th+theme_bw()+
  scale_fill_manual(values=c("grey80","cornflowerblue","firebrick4"), name="Treatment")+
  ylab("MHC I Expression Score")+xlab("")+
  theme(legend.position="none",
        axis.text = element_text(size =10),
        axis.title = element_text(size =12),
        legend.text = element_text(size =12),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))+
  geom_signif(comparisons = comp, step_increase = 0.05,tip_length = 0.01,
              size = 0.5,vjust = 0.6,
              textsize = 3,  map_signif_level = T, color="grey60")
ggsave(file=here("rna_seq/figs/","MHCI_expression_score_Treated.pdf"),width=3, height=4)




#'## R Session Info
sessionInfo()



