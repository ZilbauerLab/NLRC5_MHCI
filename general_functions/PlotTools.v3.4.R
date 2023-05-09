# PLOTTING COLOUR SCHEMES AND THEMES (shared version)
# Most recent version: fp215, 03Oct22

library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(reshape)


# IBD diagnosis colour scheme v1:

myColours_diagnosis <- c("grey75","darkblue","darkred")
colour_possibilities_diagnosis <- c("Control","CD","UC")
names(myColours_diagnosis) <- colour_possibilities_diagnosis
fillscale_diagnosis <- scale_fill_manual(name="Diagnosis", values = myColours_diagnosis, drop = T)
colscale_diagnosis <- scale_colour_manual(name="Diagnosis", values = myColours_diagnosis, drop = T)

# IBD diagnosis colour scheme v2:

myColours_diagnosis2 <- c("grey75","darkblue","darkred","purple4")
colour_possibilities_diagnosis2 <- c("Control","CD","UC","IBD-U")
names(myColours_diagnosis2) <- colour_possibilities_diagnosis2
fillscale_diagnosis2 <- scale_fill_manual(name="Diagnosis", values = myColours_diagnosis2, drop = T)
colscale_diagnosis2 <- scale_colour_manual(name="Diagnosis", values = myColours_diagnosis2, drop = T)

# IBD diagnosis colour scheme v3:

myColours_diagnosis3 <- c("grey75","purple4")
colour_possibilities_diagnosis3 <- c("Control","IBD")
names(myColours_diagnosis3) <- colour_possibilities_diagnosis3
fillscale_diagnosis3 <- scale_fill_manual(name="Diagnosis", values = myColours_diagnosis3, drop = T)
colscale_diagnosis3 <- scale_colour_manual(name="Diagnosis", values = myColours_diagnosis3, drop = T)

# IBD diagnosis colour scheme v4:

myColours_diagnosis4 <- c("grey75","darkblue")
colour_possibilities_diagnosis4 <- c("UC.Control","CD")
names(myColours_diagnosis4) <- colour_possibilities_diagnosis4
fillscale_diagnosis4 <- scale_fill_manual(name="Diagnosis", values = myColours_diagnosis4, drop = T)
colscale_diagnosis4 <- scale_colour_manual(name="Diagnosis", values = myColours_diagnosis4, drop = T)

# IBD diagnosis colour scheme v5 (RE colour scheme):

myColors_RE <- c("lightgrey","darkgoldenrod1","dodgerblue3")
color_possibilities_RE <- c("Control","UC","CD")
names(myColors_RE) <- color_possibilities_RE
fillscale_RE_diagnosis <- scale_fill_manual(name="Diagnosis", values=myColors_RE, drop=T, limits=force)
colscale_RE_diagnosis <- scale_color_manual(name="Diagnosis", values=myColors_RE, drop=T)

# IBD diagnosis colour scheme v6:

myColours_diagnosis6 <- c("grey75","darkblue","darkred","purple4","darkorange","white")
colour_possibilities_diagnosis6 <- c("Control","CD","UC","IBD-U","Other.GI","NK")
names(myColours_diagnosis6) <- colour_possibilities_diagnosis6
fillscale_diagnosis6 <- scale_fill_manual(name="Diagnosis", values = myColours_diagnosis6, drop = T)
colscale_diagnosis6 <- scale_colour_manual(name="Diagnosis", values = myColours_diagnosis6, drop = T)

# IBD diagnosis colour scheme v7 (RE colour scheme including IBDU):

myColors_RE <- c("lightgrey","darkgoldenrod1","dodgerblue3","#636B4D")
color_possibilities_RE <- c("Control","UC","CD","IBD-U")
names(myColors_RE) <- color_possibilities_RE
fillscale_RE_diagnosis2 <- scale_fill_manual(name="Diagnosis", values=myColors_RE, drop=T, limits=force)
colscale_RE_diagnosis2 <- scale_color_manual(name="Diagnosis", values=myColors_RE, drop=T)

# IBD sample site colour scheme v1:

myColours_sampsite <- c("goldenrod3","darkolivegreen3","cornflowerblue")
colour_possibilities_sampsite <- c("AC","SC","TI")
names(myColours_sampsite) <- colour_possibilities_sampsite
fillscale_sampsite <- scale_fill_manual(name="Sample Site", values = myColours_sampsite, drop = T)
colscale_sampsite <- scale_colour_manual(name="Sample Site", values = myColours_sampsite, drop = T)

# IBD sample site colour scheme v2:

myColours_sampsite2 <- c("goldenrod3","darkolivegreen3","cornflowerblue","orchid")
colour_possibilities_sampsite2 <- c("AC","SC","TI","DUO")
names(myColours_sampsite2) <- colour_possibilities_sampsite2
fillscale_sampsite2 <- scale_fill_manual(name="Sample Site", values = myColours_sampsite2, drop = T)
colscale_sampsite2 <- scale_colour_manual(name="Sample Site", values = myColours_sampsite2, drop = T)

# IBD sample site colour scheme v3:

myColours_sampsite3 <- c("goldenrod3","darkolivegreen3","cornflowerblue","orchid","darkolivegreen1","skyblue1","lightcyan")
colour_possibilities_sampsite3 <- c("AC","SC","TI","DUO","FPG","FDG","SB")
names(myColours_sampsite3) <- colour_possibilities_sampsite3
fillscale_sampsite3 <- scale_fill_manual(name="Sample Site", values = myColours_sampsite3, drop = T)
colscale_sampsite3 <- scale_colour_manual(name="Sample Site", values = myColours_sampsite3, drop = T)

# IBD inflammation colour scheme v1:

myColors_inflammation <- c("darkred","tomato","grey75","grey40")
color_possibilities_inflammation <- c("Macroscopic","Microscopic","Normal","Undetermined")
names(myColors_inflammation) <- color_possibilities_inflammation
fillscale_inflammation <- scale_fill_manual(name="Inflammation", values=myColors_inflammation, drop=T, limits=force)
colscale_inflammation <- scale_color_manual(name="Inflammation", values=myColors_inflammation, drop=T)

# IBD inflammation colour scheme v2:

myColors_inflammation2 <- c("grey75","darkred")
color_possibilities_inflammation2 <- c("N","Y")
names(myColors_inflammation2) <- color_possibilities_inflammation2
fillscale_inflammation2 <- scale_fill_manual(name="Inflammation", values=myColors_inflammation2, drop=T, limits=force)
colscale_inflammation2 <- scale_color_manual(name="Inflammation", values=myColors_inflammation2, drop=T)

# IBD sampling time point colour scheme:

myColours_timept <- c("black","purple3")
colour_possibilities_timept <- c("original","rescope")
names(myColours_timept) <- colour_possibilities_timept
fillscale_time <- scale_fill_manual(name="Time Point", values = myColours_timept, drop=T)
colscale_time <- scale_colour_manual(name="Time Point", values = myColours_timept, drop=T)

# IBD CD distribution code colour scheme:

myColours_distrib <- c("black","darkblue","goldenrod","darkolivegreen","darkred","darkorchid4","cornflowerblue")
colour_possibilities_distrib <- c("L1","L1.L4","L2","L2.L4","L3","L3.L4","L4")
names(myColours_distrib) <- colour_possibilities_distrib
fillscale_distrib <- scale_fill_manual(name="Distribution Code", values = myColours_distrib, drop=T)
colscale_distrib <- scale_colour_manual(name="Distribution Code", values = myColours_distrib, drop=T)

# IBD deep ulcer colour scheme:

myColors_ulcer <- c("grey75","black")
color_possibilities_ulcer <- c("No","Yes")
names(myColors_ulcer) <- color_possibilities_ulcer
fillscale_ulcer <- scale_fill_manual(name="Deep.Ulcer", values=myColors_ulcer, drop=T, limits=force)
colscale_ulcer <- scale_color_manual(name="Deep.Ulcer", values=myColors_ulcer, drop=T)

# IBD disease activity colour scheme:

myColors_activity <- c("black","grey50","grey75")
color_possibilities_activity <- c("active","inactive","normal")
names(myColors_activity) <- color_possibilities_activity
fillscale_activity <- scale_fill_manual(name="Disease.Activity", values=myColors_activity, drop=T, limits=force)
colscale_activity <- scale_color_manual(name="Disease.Activity", values=myColors_activity, drop=T)

# Methylation array type colour scheme:

myColours_array <- c("cyan3","magenta3")
colour_possibilities_array <- c("K450","EPIC")
names(myColours_array) <- colour_possibilities_array
fillscale_array <- scale_fill_manual(name="Array Type", values = myColours_array, drop=T)
colscale_array <- scale_colour_manual(name="Array Type", values = myColours_array, drop=T)

# Sample type colour scheme:

myColours_samptype <- c("skyblue4","yellow4")
colour_possibilities_samptype <- c("epithelium","organoid")
names(myColours_samptype) <- colour_possibilities_samptype
fillscale_samptype <- scale_fill_manual(name="Sample Type", values = myColours_samptype, drop=T)
colscale_samptype <- scale_colour_manual(name="Sample Type", values = myColours_samptype, drop=T)

# Age group colour scheme:

myColours_age <- c("goldenrod","darkorange3","firebrick")
colour_possibilities_age <- c("foetal","neonatal","paediatric")
names(myColours_age) <- colour_possibilities_age
fillscale_age <- scale_fill_manual(name="Age Group", values = myColours_age, drop = T)
colscale_age <- scale_colour_manual(name="Age Group", values = myColours_age, drop = T)

# Average MHCI clustering colour scheme:

myColors <- c("#01665e","#74AE97","#e6f5d0")
color_possibilities <- c("High_MHC","Intermediate","Low_MHC")
names(myColors) <- color_possibilities
fillscale_cluster <- scale_fill_manual(name="Cluster", values=myColors, drop=T, limits=c("High_MHC","Intermediate","Low_MHC"), na.value="white")
colscale_cluster <- scale_color_manual(name="Cluster", values=myColors, drop=T)

# AZA treatment colour scheme:

myColours <- c("palegoldenrod","darkred","goldenrod3")
colour_possibilities <- c("NT","AZA-72H","AZA-3W")
names(myColours) <- colour_possibilities
fillscale_aza <- scale_fill_manual(name="AZA Treatment", values = myColours, drop = T)
colscale_aza <- scale_colour_manual(name="AZA Treatment", values = myColours, drop = T)

# NT/TNFa/IFNg treatment colour scheme:

myColours <- c("grey75","turquoise","chartreuse4")
colour_possibilities <- c("NT","IFNg","TNFa")
names(myColours) <- colour_possibilities
fillscale_treat <- scale_fill_manual(name="Treatment", values=myColours, drop=T)
colscale_treat <- scale_colour_manual(name="Treatment", values=myColours, drop=T)

# EoE diagnosis colour scheme:

myColours_diagnosis_eoe <- c("grey75","darkred")
colour_possibilities_diagnosis_eoe <- c("Control","EoE")
names(myColours_diagnosis_eoe) <- colour_possibilities_diagnosis_eoe
fillscale_diagnosis_eoe <- scale_fill_manual(name="Diagnosis", values = myColours_diagnosis_eoe, drop = T)
colscale_diagnosis_eoe <- scale_colour_manual(name="Diagnosis", values = myColours_diagnosis_eoe, drop = T)

# EoE treatment timepoint colour scheme:

myColours_timept_eoe <- c("black","dodgerblue","dodgerblue4")
colour_possibilities_timept_eoe <- c("T0","T1","T2")
names(myColours_timept_eoe) <- colour_possibilities_timept_eoe
fillscale_time_eoe <- scale_fill_manual(name="Time Point", values=myColours_timept_eoe, drop=T)
colscale_time_eoe <- scale_colour_manual(name="Time Point", values=myColours_timept_eoe, drop=T)

# EoE diagnosis + treatment timepoint colour scheme:

myColours_eoe <- c("grey75","darkblue","darkorchid","darkred")
colour_possibilities_eoe <- c("Control","EoE.T0","EoE.T1","EoE.T2")
names(myColours_eoe) <- colour_possibilities_eoe
fillscale_eoe <- scale_fill_manual(name="DiagnosisTimepoint",values=myColours_eoe,drop=T)
colscale_eoe <- scale_colour_manual(name="DiagnosisTimepoint",values=myColours_eoe,drop=T)

# EoE diagnosis + treatment timepoint colour scheme excluding T2 timepoints (probably not necessary):

myColours_dg.tpt <- c("grey75","darkblue","darkred")
colour_possibilities_dg.tpt <- c("Control","EoE.T0","EoE.T1")
names(myColours_dg.tpt) <- colour_possibilities_dg.tpt
fillscale_dg.tpt <- scale_fill_manual(name="Diagnosis", values = myColours_dg.tpt, drop = T)
colscale_dg.tpt <- scale_colour_manual(name="Diagnosis", values = myColours_dg.tpt, drop = T)

# Heat scree plot colour scheme:

myColours_heatscree <- c("#084594","#4292c6","#9ecae1","#deebf7")
colour_possibilities_heatscree <- c("<=0.001","<=0.01","<=0.05",">0.05")
names(myColours_heatscree) <- colour_possibilities_heatscree
fillscale_heatscree <- scale_fill_manual(name="P value", values = myColours_heatscree, drop=T)
colscale_heatscree <- scale_colour_manual(name="P value", values = myColours_heatscree, drop=T)

# P-value colour scheme:

myColors_pval <- c("#084594","#4292c6","#9ecae1","#deebf7")
color_possibilities_pval<-c( "<=0.001","<=0.01","<=0.05",">0.05")
names(myColors_pval) <- color_possibilities_pval
fillscale_pval <- scale_fill_manual(name="P Value",values = myColors_pval, drop = F)


# Plot theme 1:

myTheme <- theme(axis.text=element_text(size=12), axis.title=element_text(size=14), strip.text.x = element_text(size = 12), legend.text=element_text(size=12), legend.title=element_text(size=14))

# Plot theme 2:

myTheme2 <- theme(axis.text=element_text(size=12),axis.title=element_text(size=14),strip.text.x=element_text(size=12),legend.text=element_text(size=12),legend.title=element_text(size=14))

# Plot theme 3:

ClusterTheme <- theme(legend.position="none",axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background=element_blank(),axis.line=element_blank(),panel.border=element_blank())

# Plot themes for annotation of Haberman data:

myTheme_haberman_label <- theme(legend.position="none", axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_blank(), panel.border=element_blank())
myTheme_haberman_exprs <- theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# Plot themes for average MHCI methylation plots:

myTheme_meanmeth_label <- theme(legend.position="bottom", axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_blank(), axis.line=element_blank(), panel.border=element_blank(), legend.text=element_text(size=rel(2)))
myTheme_meanmeth_meth <- theme(axis.title=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank())

# Plot theme for plotting correlation between potential covariates:

theme_correlation <- theme(axis.text=element_text(size=10, color="black"), axis.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.text=element_text(size=10), legend.title=element_text(size=10), panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())



###################
# Heat scree plot #
###################

# Heat scree plot function with categorical and continuous variables:

heat.scree <- function(Loadings, Importance, right_marg, left_marg){

    if(missing(right_marg)){
        right_marg=2.25}

    if(missing(left_marg)){
        left_marg=1}

    pca.df <- data.frame(variance=Importance, PC=seq(1:length(Importance)))
    
    scree <- ggplot(pca.df[which(pca.df$PC<=(PCs_to_view)),],aes(PC,variance))+
        geom_bar(stat="identity",color="black",fill="grey")+
        theme_bw()+
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title=element_text(size=15),plot.margin=unit(c(1.25,2,-0.2,3),"cm"))+
        ylab("Variance")+
        scale_x_continuous(breaks=seq(1,PCs_to_view,1))+
        xlab("")

# Correlate meta data with PCs and run ANOVA of each PC on each variable:
# NB, Rachel E's original script uses "spearman" to test for correlation with the continuous variables, however this can't handle ties, so updated here to "kendall"

    aov.PC.meta <- lapply(1:ncol(meta_categorical),function(covar) sapply(1:ncol(Loadings), function(PC) summary(aov(Loadings[,PC] ~ meta_categorical[,covar]))[[1]]$"Pr(>F)"[1]))

    cor.PC.meta <- lapply(1:ncol(meta_continuous), function(covar) sapply(1:ncol(Loadings), function(PC) (cor.test(Loadings[,PC], as.numeric(meta_continuous[,covar]), alternative="two.sided", method="kendall", na.action=na.omit, exact=F)$p.value)))

    names(aov.PC.meta) <- colnames(meta_categorical)
    names(cor.PC.meta) <- colnames(meta_continuous)
    aov.PC.meta <- do.call(rbind, aov.PC.meta)
    cor.PC.meta <- do.call(rbind, cor.PC.meta)
    aov.PC.meta <- rbind(aov.PC.meta, cor.PC.meta)
    aov.PC.meta <- as.data.frame(aov.PC.meta)

# Adjust:

    aov.PC.meta.adjust <- aov.PC.meta[,1:ncol(aov.PC.meta)]

# Reshape:

    avo <- aov.PC.meta.adjust[,1:PCs_to_view]
    avo.heat.no <- apply(avo,2,as.numeric)
    avo.heat <- as.data.frame(avo.heat.no)
    avo.heat$meta <- rownames(avo)
    avo.heat.melt <- melt(avo.heat, id=c("meta"))

# Cluster meta data:

    meta.var.order <- unique(avo.heat.melt$meta)[rev(ord)]
    avo.heat.melt$meta <- factor(avo.heat.melt$meta, levels=meta.var.order)

# Colour if significant:

    avo.heat.melt$Pvalue <- sapply(1:nrow(avo.heat.melt), function(x) if(avo.heat.melt$value[x]<=0.001){"<=0.001"}else{
        if(avo.heat.melt$value[x]<=0.01){"<=0.01"}else{
            if(avo.heat.melt$value[x]<=0.05){"<=0.05"}else{">0.05"}}})

    levels(avo.heat.melt$variable) <- sapply(1:PCs_to_view, function(x) paste("PC",x,sep=""))

# Create plot:

    heat <- ggplot(avo.heat.melt, aes(variable, meta, fill=Pvalue))+
        geom_tile(color="black", size=0.5)+
        theme_gray(8)+fillscale_heatscree+
        theme(axis.text=element_text(size=10,color="black"),axis.text.x=element_text(),axis.title=element_text(size=15),legend.text=element_text(size=10),legend.title=element_text(size=10),legend.position=c(1.15,0.75),legend.justification=c(1,1),plot.margin=unit(c(-0.3,right_marg,1,left_marg),"cm"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
        xlab("Principal Component")+
        ylab(NULL)

    grid.arrange(scree, heat, ncol=1,heights = c(3, 4))
}

# Heat scree plot with only categorical variables:

heat.scree.cat <- function(Loadings, Importance, right_marg, left_marg){

    if(missing(right_marg)){
        right_marg=2.25}

    if(missing(left_marg)){
        left_marg=1}

    pca.df <- data.frame(variance=Importance, PC=seq(1:length(Importance)))

    scree <- ggplot(pca.df[which(pca.df$PC<=(PCs_to_view)),],aes(PC,variance))+
        geom_bar(stat="identity",color="black",fill="grey")+
        theme_bw()+
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),axis.title=element_text(size=15),plot.margin=unit(c(1.25,2,-0.2,3),"cm"))+
        ylab("Variance")+
        scale_x_continuous(breaks=seq(1,PCs_to_view,1))+
        xlab("")

# Correlate meta data with PCs and run ANOVA of each PC on each variable:

    aov.PC.meta <- lapply(1:ncol(meta_categorical),function(covar) sapply(1:ncol(Loadings), function(PC) summary(aov(Loadings[,PC] ~ meta_categorical[,covar]))[[1]]$"Pr(>F)"[1]))
    names(aov.PC.meta) <- colnames(meta_categorical)
    aov.PC.meta <- do.call(rbind, aov.PC.meta)
    aov.PC.meta <- as.data.frame(aov.PC.meta)

# Adjust:

    aov.PC.meta.adjust <- aov.PC.meta[,1:ncol(aov.PC.meta)]

# Reshape:

    avo <- aov.PC.meta.adjust[,1:PCs_to_view]
    avo.heat.no <- apply(avo,2,as.numeric)
    avo.heat <- as.data.frame(avo.heat.no)
    avo.heat$meta <- rownames(avo)
    avo.heat.melt <- melt(avo.heat, id=c("meta"))

# Cluster meta data:

    meta.var.order <- unique(avo.heat.melt$meta)[rev(ord)]
    avo.heat.melt$meta <- factor(avo.heat.melt$meta, levels=meta.var.order)

# Colour if significant:

    avo.heat.melt$Pvalue <- sapply(1:nrow(avo.heat.melt), function(x) if(avo.heat.melt$value[x]<=0.001){"<=0.001"}else{
        if(avo.heat.melt$value[x]<=0.01){"<=0.01"}else{
            if(avo.heat.melt$value[x]<=0.05){"<=0.05"}else{">0.05"}}})

    levels(avo.heat.melt$variable) <- sapply(1:PCs_to_view, function(x) paste("PC",x,sep=""))

# Create plot:

    heat <- ggplot(avo.heat.melt, aes(variable, meta, fill=Pvalue))+
        geom_tile(color="black", size=0.5)+
        theme_gray(8)+fillscale_heatscree+
        theme(axis.text=element_text(size=10,color="black"),axis.text.x=element_text(),axis.title=element_text(size=15),legend.text=element_text(size=10),legend.title=element_text(size=10),legend.position=c(1.15,0.75),legend.justification=c(1,1),plot.margin=unit(c(-0.3,right_marg,1,left_marg),"cm"),panel.background=element_blank(),panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
        xlab("Principal Component")+
        ylab(NULL)

    grid.arrange(scree, heat, ncol=1,heights = c(3, 4))
}


################
# VOLCANO PLOT #
################

# Function created by Rachel Edgar

library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)

makeVolcano <- function(pvalue, deltabeta, dB_threshold, pval_threshold, legend_title, xlimit, ylimit){

  # Select p-values & delta beta values:

  volcano <- data.frame(Pvalue=pvalue, Delta_Beta=deltabeta)
  
  # Thresholds (delta beta and p-value)

  dB <- dB_threshold
  Pv <- pval_threshold
  
  sta_delbeta <- deltabeta[which(pvalue<=pval_threshold)] 
  sta_delbeta <- sta_delbeta[abs(sta_delbeta)>=dB]
  
  print(paste("Hypermethylated", length(sta_delbeta[which(sta_delbeta>=dB)]), sep=": "))
  print(paste("Hypomethylated", length(sta_delbeta[which(sta_delbeta<=(-dB))]) , sep=": "))
  
  volcano <- volcano[complete.cases(volcano),]
  
  # Positive delta beta is hypomethylated

  color3 <- sapply(1:nrow(volcano), function(x) if(volcano$Pvalue[x]<=Pv){
    if(abs(volcano$Delta_Beta[x])>dB){
      if(volcano$Delta_Beta[x]>dB){"Increased Methylation\n(with Potential Biological Impact)"}else{"Decreased Methylation\n (with Potential Biological Impact)"}
    }else{if(volcano$Delta_Beta[x]>0){"Increased Methylation"}else{"Decreased Methylation"}}}else{"Not Significantly Different"})
  
  volcano$Interesting_CpG3 <- color3
  
  # Define colours:

  myColors <- c(muted("red", l=80, c=30),"red",muted("blue", l=70, c=40),"blue", "grey")
  color_possibilities <-c("Decreased Methylation",
                         "Decreased Methylation\n (with Potential Biological Impact)",
                         "Increased Methylation",
                         "Increased Methylation\n(with Potential Biological Impact)",
                         "Not Significantly Different")
  names(myColors) <- color_possibilities
  colscale <- scale_color_manual(name = legend_title,
                                 values = myColors, drop = FALSE)

  # Define plot theme:
  
    th <- theme(axis.text=element_text(size=10), axis.title=element_text(size=12), strip.text = element_text(size = 12), legend.text=element_text(size=12), legend.title=element_text(size=14))

  # VOLCANO PLOT:
  
  volcano_plot <- ggplot(volcano, aes(Delta_Beta, -log10(Pvalue), color=Interesting_CpG3))+geom_point(shape=19, size=1)+theme_bw()+colscale+geom_vline(xintercept=c(-dB,dB), color="grey60")+geom_hline(yintercept=-log10(Pv), color="grey60")+ylab("P Value (-log10)")+xlab("Delta Beta")+xlim(-xlimit, xlimit)+ylim(0,ylimit)+theme(plot.margin=unit(c(1,1,1,2),"cm"))+th+guides(color = guide_legend(override.aes = list(size = 4)))
  
  # P-VALUE DISTRIBUTION PLOT:

  pval_dis <- ggplot()+geom_histogram(aes(pvalue),fill="grey90", color="black",binwidth=0.025)+theme_bw()+xlab("Nominal P Value")+th+theme(plot.margin=unit(c(1,8.8,-0.4,1.2),"cm"))+ylab("CpG Count")
  
  grid.arrange(pval_dis, volcano_plot, ncol=1,heights = c(2, 6))
  
}

