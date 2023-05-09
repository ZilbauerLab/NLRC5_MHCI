#' Functions to calculate the average expression of one or more sets of genes from bulk RNASeq or expression array data
#' Based on the Seurat AddModuleScore function: https://www.rdocumentation.org/packages/Seurat/versions/2.3.4/topics/AddModuleScore
#' 26Apr22, fp215
#'
#' @param assay.data Matrix of normalised counts or expression array data
#' @param features A list of vectors of features (each vector containing a gene or probe set of interest)
#' @param nbin Number of bins of aggregate expression levels required for all analysed features (default = 50)
#' @param controls Number of control features to be selected from the same bin per analysed feature (default = 100)
#' @param name Name preceding each average expression score produced for each vector of features
#' @seed Set a random seed (default = 1)
#' @genekey A dataframe containing a gene ID : gene symbol or probe : gene key (default = NULL)
#' @info Any additional information object to be included with output (default = NULL)
#'
#' @return Returns a list containing a vector of the input features, a vector of the control genes used, a dataframe containing the average expression measurement(s) and the additional information object if included
#' 
#' @importFrom ggplot2 cut_number
#' @importFrom Matrix rowMeans colMeans
#'

AverageExpression <- function(assay.data, features, nbin=50, controls=100, name, seed=1, genekey, info=NULL){

    cluster.length <- length(features)
    
    print(paste0("Calculating average expression for ",cluster.length," feature set(s)"))
    
    # Split dataset into nbin equal intervals based on average expression across all samples:
    
    pool <- rownames(assay.data)
    data.avg <- Matrix::rowMeans(x=assay.data[pool,])
    data.avg <- data.avg[order(data.avg)]
    data.cut <- cut_number(data.avg + rnorm(length(data.avg))/1e30, n=nbin, lables=F, right=F)
    names(data.cut) <- names(data.avg)
    
    # Calculate control gene set(s):
    
    ctrl.use <- vector(mode="list", length=cluster.length)
    
    for (i in 1:cluster.length){
        features.use <- features[[i]]
        for (j in 1:length(features.use)){
            ctrl.use[[i]] <- c(ctrl.use[[i]], names(sample(data.cut[which(data.cut==data.cut[features.use[j]])], size=controls, replace=F)))
        }
    }
    
    print(paste0("Number of control features: ",length(ctrl.use[[1]])))
    
    ctrl.scores <- matrix(data=numeric(length=1L), nrow=length(ctrl.use), ncol=ncol(assay.data))
    
    for (i in 1:length(ctrl.use)){
        features.use <- ctrl.use[[i]]
        ctrl.scores[i,] <- Matrix::colMeans(assay.data[features.use,])
    }
    
    # Calculate the average expression score:
    
    features.scores <- matrix(data=numeric(length=1L), nrow=cluster.length, ncol=ncol(assay.data))
    
    for (i in 1:cluster.length){
        features.use <- features[[i]]
        data.use <- assay.data[features.use, , drop=F]
        features.scores[i,] <- Matrix::colMeans(data.use)
    }
    
    average.expression <- features.scores - ctrl.scores
    average.expression <- as.data.frame(t(average.expression))
    rownames(average.expression) <- colnames(assay.data)
    colnames(average.expression) <- paste0(name, 1:cluster.length)
    
    # Pull together useful output:
    
    readme <- data.frame(dataset=c("FeatureSet","Controls","AverageExpression","FeatureKey","Info"),description=c("Feature set(s) of interest","List of 
features used 
to control for average feature expression","Average expression for feature set(s) of interest","Feature : gene symbol key ","Additional information provided"))
    output <- list(features,ctrl.use,average.expression,genekey,info,readme)
    names(output) <- c("FeatureSet","Controls","AverageExpression","FeatureKey","Info","README")
    return(output)
}
