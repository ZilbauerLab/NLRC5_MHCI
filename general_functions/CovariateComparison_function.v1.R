# FUNCTION TO COMPARE COVARIATES
# Input needed:
#	- meta_categorical: dataframe of categorical covariates
#	- meta_continuous: dataframe of continuous covariates
# Based on Rachel Edgar's code
# 30Sep22

library(reshape2)
library(ggplot2)


# Themes & schemes:

theme_correlation <- theme(axis.text=element_text(size=10, color="black"), axis.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.text=element_text(size=10), legend.title=element_text(size=10), panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

myColors_pval <- c("#084594","#4292c6","#9ecae1","#deebf7")
color_possibilities_pval<-c( "<=0.001","<=0.01","<=0.05",">0.05")
names(myColors_pval) <- color_possibilities_pval
fillscale_pval <- scale_fill_manual(name="P Value",values = myColors_pval, drop = F)


CompareCovars <- function(meta_categorical, meta_continuous){
    print(paste0("Comparing ", ncol(meta_categorical)+ncol(meta_continuous), " covariates"))

# Calculate correlations:

    cor_meta <- lapply(1:ncol(meta_continuous), function(con1) {
      sapply(1:ncol(meta_continuous), function(con2){
        cor.test(meta_continuous[,con1] , meta_continuous[,con2])$p.value})
    })
    cor_meta <- as.data.frame(do.call(rbind, cor_meta))
    colnames(cor_meta) <- colnames(meta_continuous)

    aov_meta <- lapply(1:ncol(meta_categorical), function(cat) {
      sapply(1:ncol(meta_continuous), function(con1){
        summary(aov(meta_continuous[, con1] ~ meta_categorical[, cat]))[[1]]$"Pr(>F)"[1] })
    })
    aov_meta <- as.data.frame(do.call(rbind, aov_meta))
    colnames(aov_meta) <- colnames(meta_continuous)

    suppressWarnings(chi_meta <- lapply(1:ncol(meta_categorical), function(cat1) {
      sapply(1:ncol(meta_categorical), function(cat2){
        chisq.test(table(meta_categorical[,cat1], meta_categorical[,cat2]))$p.value})
    }))

    chi_meta <- as.data.frame(do.call(rbind, chi_meta))
    colnames(chi_meta) <- colnames(meta_categorical)

    chi_meta <- cbind(chi_meta,aov_meta)
    con_t <- t(aov_meta)
    colnames(con_t) <- colnames(meta_categorical)
    chi_meta<-rbind(chi_meta,cbind(con_t,cor_meta))

    chi_meta$id <- c(colnames(meta_categorical),colnames(meta_continuous))
    chi_meta <- melt(chi_meta)

    chi_meta$Pvalue <- sapply(1:nrow(chi_meta), function(x) if(chi_meta$value[x]<=0.001){"<=0.001"}else{
      if(chi_meta$value[x]<=0.01){"<=0.01"}else{
        if(chi_meta$value[x]<=0.05){"<=0.05"}else{">0.05"}}})
    chi_meta$variable<-as.character(chi_meta$variable)

# Create plot input:

    p1 <- ggplot(chi_meta, aes(id, variable, fill=Pvalue))+geom_tile(color="black", size=0.5)+geom_text(aes(label=round(value, 3)))+theme_gray(8)+theme_correlation+fillscale_pval

# Output:

    return(list(covar.correlation=chi_meta,plot.input=p1))

    }
