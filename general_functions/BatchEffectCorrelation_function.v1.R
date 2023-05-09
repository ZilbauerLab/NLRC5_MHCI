# FUNCTION TO CHECK FOR CORRELATION BETWEEN META DATA AND SVA SURROGATE VARIABLES
# Input needed:
#	- meta_categorical: dataframe of categorical covariates
#	- meta_continuous: dataframe of continuous covariates
#	- meta_sva: dataframe of calculated sva surrogate variables
# 30Sep22

library(reshape2)
library(ggplot2)


# Themes & schemes:

theme_correlation <- theme(axis.text=element_text(size=10, color="black"), axis.title=element_blank(), axis.text.x=element_text(angle=90, hjust=1, vjust=0.5), legend.text=element_text(size=10), legend.title=element_text(size=10), panel.background=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())

myColors_pval <- c("#084594","#4292c6","#9ecae1","#deebf7")
color_possibilities_pval<-c( "<=0.001","<=0.01","<=0.05",">0.05")
names(myColors_pval) <- color_possibilities_pval
fillscale_pval <- scale_fill_manual(name="P Value",values = myColors_pval, drop = F)


SVcor <- function(meta_categorical, meta_continuous, meta_sva){
    print(paste0("Comparing ", ncol(meta_categorical)+ncol(meta_continuous), " variables with ", ncol(meta_sva), " sva surrogate variables"))

# Calculate correlations:

    cor_meta <- lapply(1:ncol(meta_continuous), function(con) {
      sapply(1:ncol(meta_sva), function(svs){
        cor.test(meta_continuous[,con] , meta_sva[,svs])$p.value})
    })
    cor_meta <- as.data.frame(do.call(rbind, cor_meta))
    colnames(cor_meta) <- colnames(meta_sva)

    aov_meta <- lapply(1:ncol(meta_categorical), function(cat) {
      sapply(1:ncol(meta_sva), function(svs){
        summary(aov(meta_sva[,svs] ~ meta_categorical[,cat]))[[1]]$"Pr(>F)"[1] })
    })
    aov_meta <- as.data.frame(do.call(rbind, aov_meta))
    colnames(aov_meta) <- colnames(meta_sva)

    res <- rbind(cor_meta,aov_meta)
    res$id <- c(colnames(meta_continuous),colnames(meta_categorical))
    res <- melt(res)

    res$Pvalue <- sapply(1:nrow(res), function(x) if(res$value[x]<=0.001){"<=0.001"}else{
      if(res$value[x]<=0.01){"<=0.01"}else{
        if(res$value[x]<=0.05){"<=0.05"}else{">0.05"}}})
    res$variable<-as.character(res$variable)

# Create plot input:

    p1 <- ggplot(res, aes(variable, id, fill=Pvalue))+geom_tile(color="black", size=0.5)+geom_text(aes(label=round(value, 3)))+theme_gray(8)+theme_correlation+fillscale_pval

# Output:

    return(list(sv.correlation=res,plot.input=p1))

    }
