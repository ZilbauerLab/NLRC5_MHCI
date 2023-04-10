#' ### Load Libraries
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
#library(limma)
library(here)
library(scales)
library(cowplot)
library(dplyr)
library(ggsignif)



options(stringsAsFactors = FALSE)


#' functions
source(here("general_functions","00_pretty_plots.R"))


#' ### Organoid
load(file=paste(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/"),"threebatch_combined_organoids_combatted.Rdata", sep=""))

## fix 202 and 203
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605870061_R03C01")]<-"hold"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="205605880096_R04C01")]<-"205605870061_R03C01"
epic.organoid_combined$array.id[which(epic.organoid_combined$array.id=="hold")]<-"205605880096_R04C01"

### load update diagnosis
epic.samples<-read.table(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/","AllEPICOrganoid_MHCIannotated_UpdatedSampleInfo_31Oct22.txt"), header=T, sep="\t")
epic.organoid_combined_update<-merge(epic.organoid_combined[,c("array.id" , "sample.site" ,"array.type","sentrix.pos" , 
                                                               "sample_ID","sentrix_ID","det_pval","passage.or.rescope.no_numeric", "Sample_Name",
                                                               "Wnt.type","condition","Biobank.Rachel.Replicates","batch")],
                                     epic.samples[,c("case.no", "array.id","diagnosis", "sex", "age","sample.id","numeric.passage",
                                                     "Disease.Severity","DUO.Inflammation", "TI.Inflammation", "SC.Inflammation", 
                                                     "Biologics", "Surgery","AZA", "Perianal.Disease","Avg.MHCI.CpG")],
                                     by="array.id")




#' ## delta beta correlation plot
manhattan_data<-function(CpG, pvalue, fdr, db, diagnosis, tissue, sample.type){
  data.frame(CpG=CpG, db=db, p.value=pvalue, fdr=fdr, diagnosis=diagnosis, tissue=tissue, sample.type=sample.type)
}

load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/organoid_diff_DNAm_grouped_casectrl_ThreeBatch.RData"))
diff_meth_TI_organoid<-diff_meth_TI_grouped
diff_meth_SC_organoid<-diff_meth_SC_grouped
manhattan_TI_UC_organoid<-manhattan_data(diff_meth_TI_organoid$CpG, diff_meth_TI_organoid$p.value_UC, diff_meth_TI_organoid$adjusted_p_UC, diff_meth_TI_organoid$db_ctrl_UC, "UC", "TI", "organoid")
manhattan_SC_UC_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_UC, diff_meth_SC_organoid$adjusted_p_UC, diff_meth_SC_organoid$db_ctrl_UC, "UC", "SC", "organoid")
manhattan_SC_CD_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_CD, diff_meth_SC_organoid$adjusted_p_CD, diff_meth_SC_organoid$db_ctrl_CD, "CD", "SC", "organoid")

manhattan_TI_UCCD_organoid<-manhattan_data(diff_meth_TI_organoid$CpG, diff_meth_TI_organoid$p.value_CDUC, diff_meth_TI_organoid$adjusted_p_CDUC, diff_meth_TI_organoid$db_UC_CD, "CDUC", "TI", "organoid")
manhattan_SC_UCCD_organoid<-manhattan_data(diff_meth_SC_organoid$CpG, diff_meth_SC_organoid$p.value_CDUC, diff_meth_SC_organoid$adjusted_p_CDUC, diff_meth_SC_organoid$db_UC_CD, "CDUC", "SC", "organoid")


load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data/diff_DNAm_TI_UCCTRL_organoid_ThreeBatch.RData"))
diff_meth_TI_UCCTRL_organoid$adjusted_p<-p.adjust(diff_meth_TI_UCCTRL_organoid$p.value, method = "fdr", n = nrow(diff_meth_TI_UCCTRL_organoid))
manhattan_TI_CD_organoid<-manhattan_data(rownames(diff_meth_TI_UCCTRL_organoid), diff_meth_TI_UCCTRL_organoid$p.value, diff_meth_TI_UCCTRL_organoid$adjusted_p, diff_meth_TI_UCCTRL_organoid$delta_beta, "CD", "TI", "organoid")

# #In SC compare IBD vs ctrl
# load(here("/media/redgar/Seagate Portable Drive/EBI_backup/thinkpad_backup/MHCI/IBD_bulk_integration/DNAm/data","organoid_diff_DNAm_simple_casectrl_ThreeBatch.RData"))
# manhattan_SC_UCCD_organoid<-manhattan_data(rownames(diff_meth_SC_organoid), diff_meth_SC_organoid$p.value, diff_meth_SC_organoid$adjusted_p, diff_meth_SC_organoid$delta_beta, "IBD", "SC", "organoid")




all_stats<-rbind(manhattan_TI_UC_organoid, manhattan_SC_UC_organoid, manhattan_SC_CD_organoid, 
                 manhattan_TI_CD_organoid, manhattan_TI_UCCD_organoid, manhattan_SC_UCCD_organoid)



meta<-epic.organoid_combined_update[,c("array.id","sample.site","diagnosis")]
meta$sample.site<-as.factor(meta$sample.site)
levels(meta$sample.site)<-c("DUO", "DUO", "SC" , "TI" )
meta$sample.site<-as.character(meta$sample.site)

anno_EPIC<-read.csv(here("/media/redgar/Seagate Portable Drive/EBI_backup/desktop_june2022/Documents/ibd/data","MethylationEPIC_v-1-0_B4.csv"), skip=7)



#' PLOT
CpG_plot<-function(CpG_OI, gene, segment){
  all_stats_OI<-all_stats[which(all_stats$CpG%in%CpG_OI),]
  sig<-all_stats_OI[which(all_stats_OI$fdr<0.05 & abs(all_stats_OI$db)>=0.05),]
  sig<-merge(sig, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
  
  betas_org<-as.data.frame(combat_organoid_Beta[which(rownames(combat_organoid_Beta)%in%CpG_OI),])
  betas_org$CpG<-rownames(betas_org)
  betas_org<-melt(betas_org)
  betas_org$sample.type<-"organoid"
  betas_org$sample.n<-"full"
  betas<-rbind(betas_org)
  
  betas<-merge(betas, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
  plt<-merge(betas, meta, by.x="variable", by.y="array.id")
  
  plt<-plt[which(plt$sample.site==segment),]
  
  sum_stat<-melt(tapply(plt$value, list(plt$MAPINFO,plt$diagnosis,plt$sample.site, plt$sample.type), mean))
  colnames(sum_stat)<-c("MAPINFO","diagnosis","sample.site","sample.type","value")
  sum_stat<-sum_stat[which(!(is.na(sum_stat$value))),]
  colnames(sig)[6]<-"sample.site"
  
  plt$new_color<-paste(plt$diagnosis, plt$sample.n)
  plt<-plt[which(plt$CpG%in%sig$CpG),]
  plt$diagnosis<-factor(plt$diagnosis, levels=c("Control","UC","CD"))
  
  p<-ggplot()+geom_boxplot(aes(diagnosis, value,fill=new_color, color=sample.n), plt, outlier.shape = NA)+facet_wrap(~CpG)+
    geom_jitter(aes(diagnosis, value,fill=new_color, color=sample.n), plt, width=0.25, shape=21, size=2)+
    scale_color_manual(values=c("black","grey40"))+fillscale_diagnosis_sample+theme_bw()+th+ylim(0,1.25)+ylab("DNAm Beta Value")+xlab("Diagnosis")+
    theme(legend.position="none",
          axis.title = element_text(size =12),
          legend.text = element_text(size =12),
          legend.title = element_text(size =12),
          strip.text.x = element_text(size = 12),
          strip.background = element_rect(fill="white"))+ scale_y_continuous(breaks = seq(0, 1, by = 0.25))
  
  sig_plt<-plt[,c(2,4,7,8)][!duplicated(plt[,c(2,4,7,8)]),]
  sig$sig<-"*"
  #* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001, **** FDR < 0.0001
  sig$sig<-sapply(1:nrow(sig), function(x){
    if(sig$fdr[x]<0.0001){"****"}else{
      if(sig$fdr[x]<0.001){"***"}else{
        if(sig$fdr[x]<0.01){"**"}else{
          if(sig$fdr[x]<0.05){"*"}
        }
      }
    }
  })
  
  sig_plt<-merge(sig_plt, sig[,c(1,5:7,9)],by=c("CpG","sample.type","diagnosis","sample.site"),all.x=T)
  sig_plt$sig[which(is.na(sig_plt$sig))]<-"ns"
  
  sig_plt$height<-as.factor(sig_plt$diagnosis)
  levels(sig_plt$height)<-c(1.09,1.15,1.005)
  sig_plt$height<-as.numeric(as.character(sig_plt$height))
  
  sig_plt$posistion<-as.factor(sig_plt$diagnosis)
  levels(sig_plt$posistion)<-c(2.5,1.5,2)
  sig_plt$posistion<-as.numeric(as.character(sig_plt$posistion))
  
  
  save_p<-p+
    annotate("segment", x = 1, xend = 3, y = 1, yend = 1, colour = "black") +
    annotate("segment", x = 1, xend = 1, y = 1, yend = 0.99, colour = "black") +
    annotate("segment", x = 3, xend = 3, y = 1, yend = 0.99, colour = "black") +
    annotate("segment", x = 2, xend = 3, y = 1.05, yend = 1.05, colour = "black") +
    annotate("segment", x = 2, xend = 2, y = 1.05, yend = 1.04, colour = "black") +
    annotate("segment", x = 3, xend = 3, y = 1.05, yend = 1.04, colour = "black") +
    annotate("segment", x = 1, xend = 2, y = 1.1, yend = 1.1, colour = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.1, yend = 1.09, colour = "black") +
    annotate("segment", x = 2, xend = 2, y = 1.1, yend = 1.09, colour = "black") +
    geom_text(aes(x=posistion,label = sig, y=height),data=sig_plt, size=4)
  save_p

  ggsave(file=paste(here("DNAm/figs/"),gene,"_",segment,"_DNAm_box_threebatch_updated.pdf", sep=""),save_p, width=4, height=3.5)
  ggsave(file=paste(here("DNAm/figs/jpeg/"),gene,"_",segment,"_DNAm_box_TI_threebatch_updated.jpeg", sep=""),save_p,  width=4, height=3.5)}


CpG_plot(c("cg07839457","cg07862320"), "NLRC5", "TI")
CpG_plot(c("cg07839457","cg07862320"), "NLRC5", "SC")
CpG_plot(c("cg07839457","cg07862320"), "NLRC5", "DUO")

CpG_plot(c("cg11706729","cg02756056"), "TAP1", "TI")
CpG_plot(c("cg11706729","cg02756056"), "TAP1", "SC")
CpG_plot(c("cg11706729","cg02756056"), "TAP1", "DUO")

CpG_plot(c("cg23235965","cg11594821"), "HLA_E", "TI")
CpG_plot(c("cg23235965","cg11594821"), "HLA_E", "SC")
CpG_plot(c("cg23235965","cg11594821"), "HLA_E", "DUO")

CpG_plot(c("cg23533285","cg06422189"), "HLA_B", "TI")
CpG_plot(c("cg23533285","cg06422189"), "HLA_B", "SC")
CpG_plot(c("cg23533285","cg06422189"), "HLA_B", "DUO")

CpG_plot(c("cg27537252","cg05475649"), "B2M", "TI")
CpG_plot(c("cg27537252","cg05475649"), "B2M", "SC")
CpG_plot(c("cg27537252","cg05475649"), "B2M", "DUO")

CpG_plot(c("cg18243910","cg27618777"), "TAP2", "TI")
CpG_plot(c("cg18243910","cg27618777"), "TAP2", "SC")
CpG_plot(c("cg18243910","cg27618777"), "TAP2", "DUO")




#CpG_plot(c("cg18243910","cg27618777"), "HLA_A", "TI")
#CpG_plot(c("cg18243910","cg27618777"), "HLA_A", "SC")

CpG_plot(c("cg01064373","cg18511546"), "HLA_C", "TI")
CpG_plot(c("cg01064373","cg18511546"), "HLA_C", "SC")

CpG_plot(c("cg01405582","cg11587584"), "HLA_F", "TI")
CpG_plot(c("cg01405582","cg11587584"), "HLA_F", "SC")

# CpG_plot(c("cg18243910","cg27618777"), "HLA_G", "TI")
# CpG_plot(c("cg18243910","cg27618777"), "HLA_G", "SC")

CpG_plot(c("cg24898914","cg16890093"), "PSMB8", "TI")
CpG_plot(c("cg24898914","cg16890093"), "PSMB8", "SC")

CpG_plot(c("cg06791592","cg10817441"), "PSMB9", "TI")
CpG_plot(c("cg06791592","cg10817441"), "PSMB9", "SC")

CpG_plot(c("cg15375424","cg11666365"), "IRF1", "TI")
CpG_plot(c("cg15375424","cg11666365"), "IRF1", "SC")


###############
## saving editing time
###############
CpG_plot<-function(CpG_OI, gene, segment){
  all_stats_OI<-all_stats[which(all_stats$CpG%in%CpG_OI),]
  sig<-all_stats_OI[which(all_stats_OI$fdr<0.05 & abs(all_stats_OI$db)>=0.05),]
  sig<-merge(sig, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
  
  betas_org<-as.data.frame(combat_organoid_Beta[which(rownames(combat_organoid_Beta)%in%CpG_OI),])
  betas_org$CpG<-rownames(betas_org)
  betas_org<-melt(betas_org)
  betas_org$sample.type<-"organoid"
  betas_org$sample.n<-"full"
  betas<-rbind(betas_org)
  
  betas<-merge(betas, anno_EPIC[,c("IlmnID","MAPINFO")], by.x="CpG", by.y="IlmnID")
  plt<-merge(betas, meta, by.x="variable", by.y="array.id")
  
  plt<-plt[which(plt$sample.site==segment),]
  
  sum_stat<-melt(tapply(plt$value, list(plt$MAPINFO,plt$diagnosis,plt$sample.site, plt$sample.type), mean))
  colnames(sum_stat)<-c("MAPINFO","diagnosis","sample.site","sample.type","value")
  sum_stat<-sum_stat[which(!(is.na(sum_stat$value))),]
  colnames(sig)[6]<-"sample.site"
  
  plt$new_color<-paste(plt$diagnosis, plt$sample.n)
  plt<-plt[which(plt$CpG%in%sig$CpG),]
  plt$diagnosis<-factor(plt$diagnosis, levels=c("Control","UC","CD"))
  
  p<-ggplot()+geom_boxplot(aes(diagnosis, value,fill=new_color, color=sample.n), plt, outlier.shape = NA)+facet_wrap(~CpG)+
    geom_jitter(aes(diagnosis, value,fill=new_color, color=sample.n), plt, width=0.25, shape=21, size=2)+
    scale_color_manual(values=c("black","grey40"))+fillscale_diagnosis_sample+theme_bw()+th+ylim(0,1.25)+ylab("DNAm Beta Value")+xlab("Diagnosis")+
    theme(legend.position="none",
          axis.title = element_text(size =12),
          legend.text = element_text(size =12),
          legend.title = element_text(size =12),
          strip.text.x = element_text(size = 12),
          strip.background = element_rect(fill="white"))+ scale_y_continuous(breaks = seq(0, 1, by = 0.25))
  
  sig_plt<-plt[,c(2,4,7,8)][!duplicated(plt[,c(2,4,7,8)]),]
  sig$sig<-"*"
  #* FDR < 0.05, ** FDR < 0.01, *** FDR < 0.001, **** FDR < 0.0001
  sig$sig<-sapply(1:nrow(sig), function(x){
    if(sig$fdr[x]<0.0001){"****"}else{
      if(sig$fdr[x]<0.001){"***"}else{
        if(sig$fdr[x]<0.01){"**"}else{
          if(sig$fdr[x]<0.05){"*"}
        }
      }
    }
  })
  
  sig_plt<-merge(sig_plt, sig[,c(1,5:7,9)],by=c("CpG","sample.type","diagnosis","sample.site"),all.x=T)
  sig_plt$sig[which(is.na(sig_plt$sig))]<-"ns"
  
  sig_plt$height<-as.factor(sig_plt$diagnosis)
  levels(sig_plt$height)<-c(1.09,1.15,1.005)
  sig_plt$height<-as.numeric(as.character(sig_plt$height))
  
  sig_plt$posistion<-as.factor(sig_plt$diagnosis)
  levels(sig_plt$posistion)<-c(2.5,1.5,2)
  sig_plt$posistion<-as.numeric(as.character(sig_plt$posistion))
  
  
  save_p<-p+
    annotate("segment", x = 1, xend = 3, y = 1, yend = 1, colour = "black") +
    annotate("segment", x = 1, xend = 1, y = 1, yend = 0.99, colour = "black") +
    annotate("segment", x = 3, xend = 3, y = 1, yend = 0.99, colour = "black") +
    annotate("segment", x = 2, xend = 3, y = 1.05, yend = 1.05, colour = "black") +
    annotate("segment", x = 2, xend = 2, y = 1.05, yend = 1.04, colour = "black") +
    annotate("segment", x = 3, xend = 3, y = 1.05, yend = 1.04, colour = "black") +
    annotate("segment", x = 1, xend = 2, y = 1.1, yend = 1.1, colour = "black") +
    annotate("segment", x = 1, xend = 1, y = 1.1, yend = 1.09, colour = "black") +
    annotate("segment", x = 2, xend = 2, y = 1.1, yend = 1.09, colour = "black") +
    geom_text(aes(x=posistion,label = sig, y=height),data=sig_plt, size=4)+ggtitle(gsub("_","-",gene))
  save_p}


save_p<-plot_grid(
  plot_grid(
    CpG_plot(c("cg23533285","cg06422189"), "HLA_B", "TI"),
    CpG_plot(c("cg01064373","cg18511546"), "HLA_C", "TI"),
    CpG_plot(c("cg01405582","cg11587584"), "HLA_F", "TI"),
    CpG_plot(c("cg15375424","cg11666365"), "IRF1", "TI"),
    CpG_plot(c("cg24898914","cg16890093"), "PSMB8", "TI"),
    CpG_plot(c("cg06791592","cg10817441"), "PSMB9", "TI"),
    CpG_plot(c("cg18243910","cg27618777"), "TAP2", "TI"),ncol=2),
  plot_grid(CpG_plot(c("cg23533285","cg06422189"), "HLA_B", "SC"),
    CpG_plot(c("cg01064373","cg18511546"), "HLA_C", "SC"),
    CpG_plot(c("cg01405582","cg11587584"), "HLA_F", "SC"),
    CpG_plot(c("cg15375424","cg11666365"), "IRF1", "SC"),
    CpG_plot(c("cg24898914","cg16890093"), "PSMB8", "SC"),
    CpG_plot(c("cg06791592","cg10817441"), "PSMB9", "SC"),
    CpG_plot(c("cg18243910","cg27618777"), "TAP2", "SC"),ncol=2),
  plot_grid(CpG_plot(c("cg23533285","cg06422189"), "HLA_B", "DUO"),
            CpG_plot(c("cg01064373","cg18511546"), "HLA_C", "DUO"),
            CpG_plot(c("cg01405582","cg11587584"), "HLA_F", "DUO"),
            CpG_plot(c("cg15375424","cg11666365"), "IRF1", "DUO"),
            CpG_plot(c("cg24898914","cg16890093"), "PSMB8", "DUO"),
            CpG_plot(c("cg06791592","cg10817441"), "PSMB9", "DUO"),
            CpG_plot(c("cg18243910","cg27618777"), "TAP2", "DUO"),ncol=2),ncol=3)

ggsave(file=here("DNAm/figs/","keyMHCI_supplement_DNAm_box_threebatch_updated.pdf"),save_p, width=18, height=13)




#################
### Score plot
#################
epic.organoid_combined_update$diagnosis<-as.factor(epic.organoid_combined_update$diagnosis)
epic.organoid_combined_update$diagnosis<-factor(epic.organoid_combined_update$diagnosis, levels=c("Control","UC","CD"))
epic.organoid_combined_update_TISC<-epic.organoid_combined_update[which(epic.organoid_combined_update$sample.site%in%c("TI","SC")),]

epic.organoid_combined_update_TISC$sample.site<-as.factor(epic.organoid_combined_update_TISC$sample.site)
epic.organoid_combined_update_TISC$sample.site<-factor(epic.organoid_combined_update_TISC$sample.site, levels=c("TI","SC"))

comp<-list(c("CD","Control"),c("UC","CD"),c("UC","Control"))
ggplot(epic.organoid_combined_update_TISC, aes(diagnosis, Avg.MHCI.CpG,fill=diagnosis))+geom_boxplot(outlier.shape = NA)+facet_wrap(~sample.site)+
    geom_jitter(width=0.25, shape=21, size=2)+
    fillscale_diagnosis+theme_bw()+th+ylab("Average MHC-I DNAm")+xlab("Diagnosis")+
    theme(legend.position="none",
          axis.title = element_text(size =12),
          legend.text = element_text(size =12),
          legend.title = element_text(size =12),
          strip.text.x = element_text(size = 12),
          strip.background = element_rect(fill="white"))+
  geom_signif(comparisons = comp, step_increase = 0.05,tip_length = 0.01,
              size = 0.5,vjust = 0.6,
              textsize = 3,  map_signif_level = T, color="grey60")
  

  ggsave(file=here("DNAm/figs/","MHCIscore_DNAm_box_threebatch.pdf"), width=4, height=3.5)
  ggsave(file=here("DNAm/figs/jpeg/","MHCIscore_DNAm_box_threebatch.jpeg"), save_p,  width=4, height=3.5)

