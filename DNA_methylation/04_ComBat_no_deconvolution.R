#'---
#'title: Combat array type but no deconvolution
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---


#' ### Load Libraries
library(genefilter)
library(matrixStats)
library(limma)
library(ggplot2)
library(gridExtra)
library(quadprog)
library(sva)
library(pamr)
library(limma)
library(here)

options(stringsAsFactors = FALSE)


#' ### Load Functions
source(here("general_functions","00_pretty_plots.R"))


#' ## Load Pediatric IEC data
load(here("../data/ibd_beta_botharrays_funnorm_filtered.RData"))
sampleinfo_EPIC <- read.table(here("../data/raw/EPIC_METHYLATION_DATA/Zerbino_EPICMethylationArray_SampleInfo_RE_mislabelled_fixed.txt"), header=T, sep="\t")
sampleinfo_450K <- read.table(here("../data/raw/K450_METHYLATION_DATA/Zerbino_K450MethylationArray_SampleInfo.txt"), header=T, sep="\t")
sampleinfo_EPIC<-sampleinfo_EPIC[which(sampleinfo_EPIC$array.id%in%colnames(ibd_beta_epic)),]
sampleinfo<-rbind(sampleinfo_450K,sampleinfo_EPIC)
sampleinfo$sample_ID<-sapply(1:nrow(sampleinfo), function(x) paste(sampleinfo$case.no[x],"_",sampleinfo$sample.site[x], sep=""))
sampleinfo$sentrix_ID<-gsub("\\_.*","",sampleinfo$array.id)
ibd_beta_450K<-ibd_beta_450K[which(rownames(ibd_beta_450K)%in%rownames(ibd_beta_epic)),]
ibd_beta_epic<-ibd_beta_epic[which(rownames(ibd_beta_epic)%in%rownames(ibd_beta_450K)),]
ibd_beta_epic<-ibd_beta_epic[match(rownames(ibd_beta_450K),rownames(ibd_beta_epic)),]
identical(rownames(ibd_beta_epic),rownames(ibd_beta_450K))
ibd_combo<-cbind(ibd_beta_450K,ibd_beta_epic)

sampleinfo<-sampleinfo[which(sampleinfo$array.id%in%colnames(ibd_combo)),]
sampleinfo<-sampleinfo[match(colnames(ibd_combo),sampleinfo$array.id),]
identical(sampleinfo$array.id, colnames(ibd_combo))


###############################
#' ## Combat array for not deconvoluted betas
###############################

# # impute 0 and 1
ibd_combo[ibd_combo==0]<-0.01
ibd_combo[ibd_combo==1]<-0.99

# impute NA
imputeMedianv3<-function(x) apply(x, 1, function(x){x[is.na(x)]<-median(x, na.rm=T); x}) #impute with row mean
ibd_combo<-t(imputeMedianv3(ibd_combo))

Mval<-function(beta) log2(beta/(1-beta))
edata = apply(ibd_combo, 1, Mval) # need mvalues for combat
edata = as.data.frame(edata)
edata = t(edata)



#mod = model.matrix(~as.factor(Tissue), data=tcell_sampleinfo) # can not protect cause confounded
batch = sampleinfo$array.type
combat_ibd_mval = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE)

#Back to betas
betas<-function(M) 2^M/((2^M)+1)
combat_ibd_Beta = apply(combat_ibd_mval, 1, betas) # need mvalues for combat
combat_ibd_Beta = as.data.frame(combat_ibd_Beta)
combat_ibd_Beta = t(combat_ibd_Beta)
combat_ibd_Beta<-as.data.frame(combat_ibd_Beta)

combat_ibd_Beta<-t(imputeMedianv3(combat_ibd_Beta))
sampleinfo_DNAm<-sampleinfo
save(combat_ibd_Beta,sampleinfo_DNAm, file=paste(here("DNAm/data/"),"ibd_not_adjusted_combatted.Rdata", sep=""))



