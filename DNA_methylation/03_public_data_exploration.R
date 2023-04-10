#'---
#'title: GEO data for comparison
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---

#' ### Load Libraries
suppressMessages(library(minfi))
suppressMessages(library(reshape))
library(ggplot2)
library(RColorBrewer)
library(here)
library(GEOquery)
library(GEOmetadb)
library(sva)
library(pamr)
library(limma)

source(here("general_functions","00_pretty_plots.R"))
options(stringsAsFactors = FALSE)



#' ### Download Fetal raw data
#wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3703/E-MTAB-3703.raw.1.zip
#wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-3703/E-MTAB-3703.sdrf.txt
path<-"DNAm/data/public/fetal_gut/E-MTAB-3703.raw.1"
sampleinfo_fetal <- read.table(here("DNAm/data/public/fetal_gut/", "E-MTAB-3703.sdrf.txt"), header=T, sep="\t")

#sample info duplicated with red and green
sampleinfo_fetal<-sampleinfo_fetal[seq(1, nrow(sampleinfo_fetal), 2),c(1:4,15)]
sampleinfo_fetal$Assay.Name<-gsub("_Grn","",sampleinfo_fetal$Assay.Name)
sampleinfo_fetal$region<-sapply(1:nrow(sampleinfo_fetal), function(x) strsplit(as.character(sampleinfo_fetal$Characteristics.cell.type.)[x], " ")[[1]][4])
sampleinfo_fetal$Characteristics.organism.<-NULL
sampleinfo_fetal$Characteristics.cell.type.<-NULL
sampleinfo_fetal$sentrix<-sapply(1:nrow(sampleinfo_fetal), function(x) strsplit(as.character(sampleinfo_fetal$Assay.Name)[x], "_")[[1]][1])
sampleinfo_fetal$sample.site<-paste(sampleinfo_fetal$Characteristics.developmental.stage.,sampleinfo_fetal$region)
sampleinfo_fetal$array.id.path <- file.path(here(path, sampleinfo_fetal$Assay.Name))

rgset_fetal <- read.metharray(sampleinfo_fetal$array.id.path, verbose = TRUE)

# Background and dye bias correction with noob thhrough funnorm implemented in minfi
# http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_fetal)
fetal_beta<-getBeta(MSet.illumina)


#' ### Load Pediatric data
load(here("../data/ibd_beta_botharrays_funnorm_filtered.RData"))

###combined meta 
sampleinfo_EPIC <- read.table(here("../data/raw/EPIC_METHYLATION_DATA/Zerbino_EPICMethylationArray_SampleInfo_RE_mislabelled_fixed.txt"), header=T, sep="\t")
sampleinfo_450K <- read.table(here("../data/raw/K450_METHYLATION_DATA/Zerbino_K450MethylationArray_SampleInfo.txt"), header=T, sep="\t")

sampleinfo_EPIC<-sampleinfo_EPIC[which(sampleinfo_EPIC$array.id%in%colnames(ibd_beta_epic)),]

sampleinfo<-rbind(sampleinfo_450K,sampleinfo_EPIC)
sampleinfo$sample_ID<-sapply(1:nrow(sampleinfo), function(x) paste(sampleinfo$case.no[x],"_",sampleinfo$sample.site[x], sep=""))
sampleinfo$sentrix_ID<-gsub("\\_.*","",sampleinfo$array.id)

### overlapping probes
ibd_beta_450K<-ibd_beta_450K[which(rownames(ibd_beta_450K)%in%rownames(ibd_beta_epic)),]
ibd_beta_epic<-ibd_beta_epic[which(rownames(ibd_beta_epic)%in%rownames(ibd_beta_450K)),]
ibd_beta_epic<-ibd_beta_epic[match(rownames(ibd_beta_450K),rownames(ibd_beta_epic)),]
identical(rownames(ibd_beta_epic),rownames(ibd_beta_450K))
ibd_combo<-cbind(ibd_beta_450K,ibd_beta_epic)

sampleinfo<-sampleinfo[which(sampleinfo$array.id%in%colnames(ibd_combo)),]
sampleinfo<-sampleinfo[match(colnames(ibd_combo),sampleinfo$array.id),]
identical(sampleinfo$array.id, colnames(ibd_combo))

### IBD vs controls in each tissue
sampleinfo$diagnosis_simple<-as.factor(sampleinfo$diagnosis)
levels(sampleinfo$diagnosis_simple)<-c("IBD","Control","IBD","IBD","IBD")

### More specific grouping
sampleinfo$diagnosis_grouped<-as.factor(sampleinfo$diagnosis)
levels(sampleinfo$diagnosis_grouped)<-c("CD","Control","CD","UC","UC")
sampleinfo$diagnosis_grouped<-factor(sampleinfo$diagnosis_grouped, levels=c("Control", "CD", "UC"))

###  More simple grouping of samplesite
sampleinfo$sample.site_grouped<-as.factor(sampleinfo$sample.site)
levels(sampleinfo$sample.site_grouped)<-c("colon","colon","ileum")     









#' ## Combine the Fetal the Pediatric
fetal_beta<-fetal_beta[which(rownames(fetal_beta)%in%rownames(ibd_combo)),]
ibd_combo_joint<-ibd_combo[which(rownames(ibd_combo)%in%rownames(fetal_beta)),]
fetal_beta<-fetal_beta[match(rownames(ibd_combo_joint),rownames(fetal_beta)),]
identical(rownames(ibd_combo_joint), rownames(fetal_beta))

ibd_fetal_beta<-cbind(ibd_combo_joint, fetal_beta)

#PCA with fetal
pca_res_joint <- prcomp(t(ibd_fetal_beta))
Loadings_joint<-as.data.frame(pca_res_joint$x)
Loadings_joint$array.id<-rownames(Loadings_joint)
Loadings_joint$tissue_plt<-c(as.character(sampleinfo$sample.site),sampleinfo_fetal$region)
Loadings_joint$stage<-c(rep("paediatric", nrow(sampleinfo)),sampleinfo_fetal$Characteristics.developmental.stage.)

ggplot(Loadings_joint, aes(PC1, PC2, fill=tissue_plt, color=stage))+geom_point(shape=21,size=2)+theme_bw()+
  scale_fill_manual(values=c("#1a9850","#006d2c","#a6d96a","dodgerblue4","cornflowerblue"))+scale_color_manual(values=c("black", "white"))

ggsave(here("DNAm/figs","PC1_PC2_fetal.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_fetal.jpeg"), width = 7.5, height = 6)









#################
#'## Pulling GEO Datasets
################


### Pull all 450K studies under GPL13534 
getSQLiteFile(destdir=here("DNAm/data/public/"))
con <- dbConnect(SQLite(), here("DNAm/data/public/GEOmetadb.sqlite"))
dbListFields(con, "gsm")
x<-dbGetQuery(con, "select title,description,series_id,gsm,source_name_ch1,characteristics_ch1 from gsm where gpl='GPL13534' OR gpl='GPL21145'")
meta<-x
print(paste("Samples: ",nrow(meta), "   Studies: ", length(unique(meta$series_id)), sep=""))


### search for terms that tissue samples would likely be under
terms<-"gut|Gut|colon|Colon|intestine|Intestine"
gut_meta<-meta[unique(c(grep(terms,meta$source_name_ch1),
                        grep(terms,meta$title), 
                        grep("terms",meta$characteristics_ch1))) ,]
print(paste("Gut Samples: ",nrow(gut_meta), "   Gut Studies: ", length(unique(gut_meta$series_id)), sep=""))


#Exclude samples that may not actually be the tissue of interest, non-human samples and cancer
#These exclusion terms can be generated by looking through the terms present for terms you don't want
# unique(gut_meta$source_name_ch1)
# unique(gut_meta$characteristics_ch1)
# unique(gut_meta$title)

gut_meta<-gut_meta[-grep("Guthrie|FFPE|adenocarcinoma|after treatment|tumor|cancer", gut_meta$source_name_ch1),]
gut_meta<-gut_meta[-grep("cancer|Tumor|HCT116", gut_meta$title),]
gut_meta<-gut_meta[-grep("GSE57342|GSE53051|GSE51815|GSE83520", gut_meta$series_id),]

print(paste("Gut Samples: ",nrow(gut_meta), "   Gut Studies: ", length(unique(gut_meta$series_id)), sep=""))

# any ibd only from colon mucosa
gut_ibd_meta<-gut_meta[grep("ibd|IBD|UC|CD|crohn|colitis|Crohn|Colitis", gut_meta$source_name_ch1),]
gut_ibd_meta<-gut_meta[grep("ibd|IBD|UC|CD|crohn|colitis|Crohn|Colitis", gut_meta$characteristics_ch1),]
gut_ibd_meta<-gut_meta[grep("ibd|IBD|UC|CD|crohn|colitis|Crohn|Colitis", gut_meta$title),]


# with age
unique(gut_meta[-grep("age|Age", gut_meta$characteristics_ch1),]$series_id)
unique(gut_meta[grep("age|Age", gut_meta$characteristics_ch1),]$series_id)

tapply(gut_meta$gsm, gut_meta$series_id, length)
tapply(gut_meta[grep("age|Age", gut_meta$characteristics_ch1),]$gsm, gut_meta[grep("age|Age", gut_meta$characteristics_ch1),]$series_id, length)

gut_meta$age<-sapply(1:nrow(gut_meta), function(x) {
  chara<-strsplit(gut_meta$characteristics_ch1[x], "\t|;")[[1]]
  gsub(".*: ","",chara[grep("age",chara)][1])
})

fetal<-gut_meta[grep("fetal",gut_meta$source_name_ch1),]

# fetal geo
GSE116754<- as.data.frame(exprs(getGEO("GSE116754")[[1]]))
GSE116754_fetal_gut<-GSE116754[,which(colnames(GSE116754)%in%fetal$gsm)]
save(GSE116754_fetal_gut,fetal, file=paste(here("DNAm/data/public/fetal_gut"),"GSE116754_fetal.RData",sep=""))



#load(here("../../ibd/data/public_other/fetal_gut","GSE116754_fetal.RData")) # only in old repo

# combine with other fetal cohort and pediatric
GSE116754_fetal_gut<-GSE116754_fetal_gut[which(rownames(GSE116754_fetal_gut)%in%rownames(ibd_fetal_beta)),]
ibd_fetal_beta<-ibd_fetal_beta[which(rownames(ibd_fetal_beta)%in%rownames(GSE116754_fetal_gut)),]
GSE116754_fetal_gut<-GSE116754_fetal_gut[match(rownames(ibd_fetal_beta),rownames(GSE116754_fetal_gut)),]
identical(rownames(ibd_fetal_beta), rownames(GSE116754_fetal_gut))

ibd_combo_AE_GEO<-cbind(ibd_fetal_beta, GSE116754_fetal_gut)

pca_res_joint <- prcomp(t(ibd_combo_AE_GEO[complete.cases(ibd_combo_AE_GEO), ]))
Loadings_joint<-as.data.frame(pca_res_joint$x)
Loadings_joint$array.id<-rownames(Loadings_joint)
Loadings_joint$tissue_plt<-c(as.character(sampleinfo$sample.site),sampleinfo_fetal$region, rep("small", 3))
Loadings_joint$stage<-c(rep("paediatric", nrow(sampleinfo)),sampleinfo_fetal$Characteristics.developmental.stage., rep("foetal",3))

ggplot(Loadings_joint, aes(PC1, PC2, fill=tissue_plt, color=stage))+geom_point(shape=21,size=2)+theme_bw()+
  scale_fill_manual(values=c("#1a9850","#006d2c","#a6d96a","dodgerblue4","cornflowerblue"))+scale_color_manual(values=c("black", "white"))

ggsave(here("DNAm/figs","PC1_PC2_fetal_GEO.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_fetal_GEO.jpeg"), width = 7.5, height = 6)


## Combine just the public fetal datasets
fetal_beta<-fetal_beta[which(rownames(fetal_beta)%in%rownames(GSE116754_fetal_gut)),]
GSE116754_fetal_gut<-GSE116754_fetal_gut[which(rownames(GSE116754_fetal_gut)%in%rownames(fetal_beta)),]
fetal_beta<-fetal_beta[match(rownames(GSE116754_fetal_gut),rownames(fetal_beta)),]
identical(rownames(fetal_beta), rownames(GSE116754_fetal_gut))

fetal_public<-cbind(fetal_beta,GSE116754_fetal_gut)



#################################### 
#' ## Loading public data (but not fetal) for more PCA
#################################### 
#'## is the PCA pattern t cell contamination??

GSE50222<-meta[grep("GSE50222", meta$series_id),]
#just one time point per sample
GSE50222<-GSE50222[grep("outside",GSE50222$description),]
GSE50222_beta<- as.data.frame(exprs(getGEO("GSE50222")[[1]]))
GSE50222_beta<-GSE50222_beta[,which(colnames(GSE50222_beta)%in%GSE50222$gsm)]
dim(GSE50222_beta)
save(GSE50222_beta, file=paste(here("DNAm/data/public/"),"GSE50222_beta.RData",sep=""))


#load(here("DNAm/data/public","GSE50222_beta.RData"))

GSE50222_beta<-GSE50222_beta[which(rownames(GSE50222_beta)%in%rownames(ibd_combo)),]
ibd_combo_joint<-ibd_combo[which(rownames(ibd_combo)%in%rownames(GSE50222_beta)),]
GSE50222_beta<-GSE50222_beta[match(rownames(ibd_combo_joint),rownames(GSE50222_beta)),]
identical(rownames(ibd_combo_joint), rownames(GSE50222_beta))

GEO_ibd_beta<-cbind(ibd_combo_joint, GSE50222_beta)

#' ## Joint combined PCA
pca_res_all <- prcomp(t(ibd_combo_joint))
Loadings_all<-as.data.frame(pca_res_all$x)
Loadings_all$array.id<-rownames(Loadings_all)
Loadings_all_meta<-merge(Loadings_all, sampleinfo, by="array.id")
ggplot(Loadings_all_meta, aes(PC1, PC2, fill=diagnosis_grouped, color=sample.site_grouped))+geom_point(shape=21,size=2)+theme_bw()+fillscale_diagnosis+scale_color_manual(values=c("black","white"))


#' PCA with tcell
pca_res_joint <- prcomp(t(GEO_ibd_beta[complete.cases(GEO_ibd_beta),]))
Loadings_joint<-as.data.frame(pca_res_joint$x)
Loadings_joint$array.id<-rownames(Loadings_joint)
Loadings_joint$tissue_plt<-c(as.character(sampleinfo$sample.site),rep("T-cell", 16))

ggplot(Loadings_joint, aes(PC1, PC2, fill=tissue_plt))+geom_point(shape=21,size=2)+theme_bw()+
  scale_fill_manual(values=c("#1a9850","#a6d96a","#e31a1c","cornflowerblue"))

ggsave(here("DNAm/figs","PC1_PC2_GEO.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_GEO.jpeg"), width = 7.5, height = 6)






##################### 
#' ## more tissues
##################### 
GSE31848<- as.data.frame(exprs(getGEO("GSE31848")[[1]]))
GSE31848_meta<-meta[grep("GSE31848", meta$series_id),]
GSE31848_meta$tissue<-sapply(1:nrow(GSE31848_meta), function(x) gsub("tissue type: ","",strsplit(GSE31848_meta$characteristics_ch1[x], ";\t")[[1]][4]))
GSE31848_meta$type<-sapply(1:nrow(GSE31848_meta), function(x) gsub("cell type: ","",strsplit(GSE31848_meta$characteristics_ch1[x], ";\t")[[1]][3]))

GSE31848_meta_somatic<-GSE31848_meta[which(GSE31848_meta$tissue%in%c("Adipose","Blood","Brain","Kidney","Lung", "Stomach")),]
GSE31848_somatic<-GSE31848[,which(colnames(GSE31848)%in%GSE31848_meta_somatic$gsm)]
GSE31848_somatic<-GSE31848_somatic[which(rownames(GSE31848_somatic)%in%rownames(GEO_ibd_beta)),]
GEO_ibd_beta_multi<-GEO_ibd_beta[which(rownames(GEO_ibd_beta)%in%rownames(GSE31848_somatic)),]
GSE31848_somatic<-GSE31848_somatic[match(rownames(GEO_ibd_beta_multi),rownames(GSE31848_somatic)),]
identical(rownames(GEO_ibd_beta_multi), rownames(GSE31848_somatic))

GEO_ibd_beta_multi<-cbind(GEO_ibd_beta_multi,GSE31848_somatic)


#' ## PCA with tcell plus more
pca_res_GEO <- prcomp(t(GEO_ibd_beta_multi[complete.cases(GEO_ibd_beta_multi),]))
Loadings_GEO<-as.data.frame(pca_res_GEO$x)
Loadings_GEO$array.id<-rownames(Loadings_GEO)
Loadings_GEO$tissue_plt<-c(as.character(sampleinfo$sample.site),rep("T-cell", 16), GSE31848_meta_somatic$tissue)

ggplot(Loadings_GEO, aes(PC1, PC2, fill=tissue_plt))+geom_point(shape=21,size=2)+theme_bw()+
  scale_fill_manual(values=c("#1a9850","#d4b9da","#fd8d3c","#bcbddc","#8c96c6","#88419d","#a6d96a","#6e016b","#e31a1c","cornflowerblue"))


##################### 
#' ## even more tissues
##################### 
GSE48472_meta<-meta[grep("GSE48472", meta$series_id),]
GSE48472<- as.data.frame(exprs(getGEO("GSE48472")[[1]]))
identical(colnames(GSE48472), GSE48472_meta$gsm)

GSE48472_meta$tissue<-sapply(1:nrow(GSE48472_meta), function(x) gsub("tissue: ","",strsplit(GSE48472_meta$characteristics_ch1[x], ";\t")[[1]][2]))
GSE48472_meta$tissue[which(GSE48472_meta$tissue=="scfat")]<-"Fat"
GSE48472_meta$tissue[which(GSE48472_meta$tissue=="blood")]<-"Blood"

GSE48472<-GSE48472[which(rownames(GSE48472)%in%rownames(GEO_ibd_beta_multi)),]
GEO_ibd_beta_multi<-GEO_ibd_beta_multi[which(rownames(GEO_ibd_beta_multi)%in%rownames(GSE48472)),]
GSE48472<-GSE48472[match(rownames(GEO_ibd_beta_multi),rownames(GSE48472)),]
identical(rownames(GEO_ibd_beta_multi), rownames(GSE48472))

GEO_ibd_beta_multi<-cbind(GEO_ibd_beta_multi,GSE48472)

#' save all for combined analysis
save(GEO_ibd_beta_multi, sampleinfo,GSE31848_meta_somatic,GSE48472_meta, file=paste(here("DNAm/data/public/"), "ibd_GEO_combined.RData",sep=""))

#load(here("DNAm/data/public/","ibd_GEO_combined.RData"))

dim(GEO_ibd_beta_multi)

#PCA with tcell plus more
pca_res_GEO <- prcomp(t(GEO_ibd_beta_multi[complete.cases(GEO_ibd_beta_multi),]))
Loadings_GEO<-as.data.frame(pca_res_GEO$x)
Loadings_GEO$array.id<-rownames(Loadings_GEO)
Loadings_GEO$tissue_plt<-c(as.character(sampleinfo$sample.site),rep("T-cell", 16), GSE31848_meta_somatic$tissue, GSE48472_meta$tissue)
Loadings_GEO$study<-c(as.character(sampleinfo$diagnosis_simple),rep("Public", 99))


ggplot(Loadings_GEO, aes(PC1, PC2, fill=tissue_plt, color=study))+geom_point(shape=21,size=2)+theme_bw()+
  scale_fill_manual(values=c("#1a9850","#fa9fb5","#fd8d3c","#9ebcda","#fed976","#fa9fb5","#fed976","#8856a7","#8856a7","grey","grey","#beaed4","#8856a7","#d95f02",
                             "#a6d96a","#8856a7","#8856a7","#e31a1c","cornflowerblue"))+
  scale_color_manual(values=c("black", "grey50", "white"))

ggsave(here("DNAm/figs","PC1_PC2_botharrays_multi_tissue.pdf"), width = 7.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_botharrays_multi_tissue.jpeg"), width = 7.5, height = 6)



######### #PCA with tcell plus more,GEO_ibd_beta_multi from other script on tissue
load(here("DNAm/data/public","ibd_GEO_combined.RData"))


fetal_public<-fetal_public[which(rownames(fetal_public)%in%rownames(GEO_ibd_beta_multi)),]
GEO_ibd_beta_multi<-GEO_ibd_beta_multi[which(rownames(GEO_ibd_beta_multi)%in%rownames(fetal_public)),]
fetal_public<-fetal_public[match(rownames(GEO_ibd_beta_multi),rownames(fetal_public)),]
identical(rownames(fetal_public), rownames(GEO_ibd_beta_multi))

GEO_ibd_beta_multi_fetal<-cbind(GEO_ibd_beta_multi,fetal_public)


pca_res_GEO <- prcomp(t(GEO_ibd_beta_multi_fetal[complete.cases(GEO_ibd_beta_multi_fetal),]))
Loadings_GEO<-as.data.frame(pca_res_GEO$x)
Loadings_GEO$array.id<-rownames(Loadings_GEO)
Loadings_GEO$tissue_plt<-c(as.character(sampleinfo$sample.site),rep("T-cell", 16), GSE31848_meta_somatic$tissue, GSE48472_meta$tissue, sampleinfo_fetal$region, rep("small", 3))
Loadings_GEO$study<-c(as.character(sampleinfo$diagnosis_simple),rep("Public", 114))


GSE31848_meta_somatic$age<-sapply(1:nrow(GSE31848_meta_somatic), function(x) strsplit(GSE31848_meta_somatic$characteristics_ch1[x], "\t")[[1]][5])
GSE31848_meta_somatic$age<-gsub("fetal vs adult tissue: |;", "",GSE31848_meta_somatic$age)

Loadings_GEO$age<-c(as.character(sampleinfo$age.group),rep("Adult", 16), GSE31848_meta_somatic$age, 
                    rep("Adult", 56), sampleinfo_fetal$Characteristics.developmental.stage., rep("Fetal", 3))




ggplot(Loadings_GEO, aes(PC1, PC2, fill=tissue_plt, color=study))+geom_point(shape=21,size=2)+theme_bw()+
  scale_fill_manual(values=c("#1a9850","#fa9fb5","#fd8d3c","#9ebcda","#fed976","#fa9fb5","#fed976","#8856a7","#a6d96a","#8856a7","grey","grey","#beaed4","#8856a7","#d95f02",
                             "#a6d96a","#006d2c","#8856a7","#8856a7","#e31a1c","cornflowerblue"))+
  scale_color_manual(values=c("black", "grey50", "white"))

Loadings_GEO$age<-as.factor(Loadings_GEO$age)
levels(Loadings_GEO$age)<-c( "Adult","Fetal","Fetal","Paediatric")
ggplot(Loadings_GEO, aes(PC1, PC2, fill=tissue_plt, color=age))+geom_point(shape=21,size=2.5)+theme_bw()+
  scale_fill_manual(values=c("#1a9850","#fa9fb5","#fd8d3c","#9ebcda","#fed976","#fa9fb5","#fed976","#8856a7","#006d2c","#8856a7","grey","grey","#beaed4","#8856a7","#d95f02",
                             "#a6d96a","dodgerblue4","#8856a7","#8856a7","#e31a1c","cornflowerblue"))+
  scale_color_manual(values=c("white", "black","white"))


unique(Loadings_GEO$tissue_plt)[order(unique(Loadings_GEO$tissue_plt))]

ggplot(Loadings_GEO, aes(PC1, PC2, fill=tissue_plt, color=age))+geom_point(shape=21,size=3)+theme_bw()+
  scale_fill_manual(values=c("#1a9850","#fcbba1","#cb181d","#F6EFBF","#F6EFBF","#fcbba1","#F6EFBF",
                             "#fcbba1","#1a9850","#E1DEE8","#E1DEE8","#fcbba1","#E1DEE8","#E1DEE8","#cb181d",
                             "#a6d96a","cornflowerblue",
                             "#fcbba1","#E1DEE8","#cb181d",
                             "cornflowerblue"), name="Tissue")+
  scale_color_manual(values=c("white", "black","grey70"), name="Age Group")+th


ggsave(here("DNAm/figs","PC1_PC2_fetal_GEO_multitissue.pdf"), width = 8.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_fetal_GEO_multitissue.jpeg"), width = 8.5, height = 6)


table(Loadings_GEO$study)
table(Loadings_GEO$tissue_plt)
table(Loadings_GEO$tissue_plt,Loadings_GEO$age)



#############
#' ## organoids to include in the PCA
#############

#' wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.raw.1.zip
#' wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.raw.2.zip
#' wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.raw.3.zip
#' wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.raw.4.zip
#' wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.raw.5.zip
#' mkdir E-MTAB-4957.raw
#' unzip E-MTAB-4957.raw.1.zip -d E-MTAB-4957.raw
#' unzip E-MTAB-4957.raw.2.zip -d E-MTAB-4957.raw
#' unzip E-MTAB-4957.raw.3.zip -d E-MTAB-4957.raw
#' unzip E-MTAB-4957.raw.4.zip -d E-MTAB-4957.raw
#' unzip E-MTAB-4957.raw.5.zip -d E-MTAB-4957.raw

#' wget https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-4957/E-MTAB-4957.sdrf.txt
path<-"DNAm/data/public/organoid/E-MTAB-4957.raw"

sampleinfo_organoid <- read.table(here("DNAm/data/public/organoid/E-MTAB-4957.sdrf.txt"), header=T, sep="\t")
sampleinfo_organoid<-sampleinfo_organoid[,c(1:15,30:34,38:41)]

#sample info cleanup
sampleinfo_organoid$sentrix<-sapply(1:nrow(sampleinfo_organoid), function(x) strsplit(as.character(sampleinfo_organoid$Assay.Name)[x], "_")[[1]][1])
sampleinfo_organoid$array.id.path <- file.path(here(path, sampleinfo_organoid$Assay.Name))
rgset_organoid <- read.metharray(sampleinfo_organoid$array.id.path, verbose = TRUE)

# Background and dye bias correction with noob thhrough funnorm implemented in minfi
#http://bioconductor.org/help/course-materials/2015/BioC2015/methylation450k.html
MSet.illumina <- preprocessFunnorm(rgset_organoid)
organoid_beta<-getBeta(MSet.illumina)


#' ## add organoids to the rest of data
organoid_beta<-organoid_beta[which(rownames(organoid_beta)%in%rownames(GEO_ibd_beta_multi_fetal)),]

organoid_beta<-organoid_beta[match(rownames(GEO_ibd_beta_multi_fetal), rownames(organoid_beta)),]
identical(rownames(organoid_beta),rownames(GEO_ibd_beta_multi_fetal))
save(organoid_beta, file=paste(here("DNAm/data/public/organoid/"),"organoid_beta.Rdata", sep=""))

#' ## Combine everything 
GEO_ibd_beta_multi_fetal_organoid<-cbind(GEO_ibd_beta_multi_fetal, organoid_beta)
dim(GEO_ibd_beta_multi_fetal_organoid)

save(GEO_ibd_beta_multi_fetal_organoid, file=paste(here("DNAm/data/public/"),"GEO_ibd_beta_multi_fetal_organoid.RData", sep=""))

#' # Fully combined PCA
pca_res_GEO <- prcomp(t(GEO_ibd_beta_multi_fetal_organoid[complete.cases(GEO_ibd_beta_multi_fetal_organoid),]))
Loadings_GEO<-as.data.frame(pca_res_GEO$x)
Loadings_GEO$array.id<-rownames(Loadings_GEO)
Loadings_GEO$tissue_plt<-c(as.character(sampleinfo$sample.site),rep("T-cell", 16), GSE31848_meta_somatic$tissue, 
                           GSE48472_meta$tissue, sampleinfo_fetal$region, rep("small", 3), 
                           paste(sampleinfo_organoid$Characteristics.biosource.type., sampleinfo_organoid$Factor.Value.sampling.site.))

Loadings_GEO$study<-c(as.character(sampleinfo$diagnosis_simple),rep("Public", 248))

Loadings_GEO$age<-c(as.character(sampleinfo$age.group),rep("Adult", 16), GSE31848_meta_somatic$age, 
                    rep("Adult", 56), sampleinfo_fetal$Characteristics.developmental.stage., rep("Fetal", 3),
                    sampleinfo_organoid$Factor.Value.developmental.stage.)
              
              # save(Loadings_GEO, file="../../data/PCA_with_organoids.RData")
              # 
              # 
              # 
              # 
              # #### run on serveer load here for plotting
              # load("../../data/PCA_with_organoids.RData")

unique(Loadings_GEO$tissue_plt[order(Loadings_GEO$tissue_plt)])

Loadings_GEO$condition<-Loadings_GEO$tissue_plt
Loadings_GEO$condition[grep("organoid",Loadings_GEO$tissue_plt)]<-"Organoid"
Loadings_GEO$condition[which(Loadings_GEO$condition!="Organoid")]<-"Primary"

Loadings_GEO$age[grep("juvenile|paediatric",Loadings_GEO$age)]<-"Paediatric"
Loadings_GEO$age[grep("fetal|foetal",Loadings_GEO$age)]<-"Fetal"
Loadings_GEO$age[grep("Adult|adult",Loadings_GEO$age)]<-"Adult"


ggplot(Loadings_GEO, aes(PC1, PC2, fill=tissue_plt, color=condition, shape=age))+geom_point(size=3)+theme_bw()+
  scale_fill_manual(values=c("#1a9850","#fcbba1","#cb181d","#F6EFBF","#F6EFBF",
                             "#fcbba1","#F6EFBF","#fcbba1","#1a9850","#E1DEE8",
                             "#E1DEE8","#fcbba1","#E1DEE8", "#a6d96a","#1a9850", 
                             "green", "cornflowerblue", "#a6d96a", "#a6d96a","#E1DEE8",
                             "cornflowerblue","#E1DEE8","#1a9850","#1a9850", 
                             "cornflowerblue", "#a6d96a", "cornflowerblue","#cb181d",
                             "#a6d96a","cornflowerblue",
                             "#fcbba1","#E1DEE8","#cb181d","cornflowerblue"), name="Tissue")+
  scale_color_manual(values=c("black", "white"), name="Age Group")+th+
  scale_shape_manual(values = c(21, 22, 21))+
  guides(fill = guide_legend(override.aes=list(shape=21)),
         color = guide_legend(override.aes=list(shape=21)))

ggsave(here("DNAm/figs","PC1_PC2_GEO_organoid_fetal.pdf"), width = 14, height = 10)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_GEO_organoid_fetal.jpeg"), width = 14, height = 10)





#' ## Spread is between blood and organoids with other inbetween
Loadings_GEO_less<-Loadings_GEO[grep("T-cell|small|organoid|large|SC|Blood|purified|organism|TI|AC",Loadings_GEO$tissue_plt),]
Loadings_GEO_less<-Loadings_GEO_less[-grep("stomach",Loadings_GEO_less$tissue_plt),]

unique(Loadings_GEO_less$tissue_plt[order(Loadings_GEO_less$tissue_plt)])
Loadings_GEO_less$tissue_plt<-as.factor(Loadings_GEO_less$tissue_plt)

levels(Loadings_GEO_less$tissue_plt)<-c("Colon","Whole Blood", "Colon","Organoid - Colon","Organoid - Colon","Organoid - Other GI",
                                        "Organoid - Small Intestine", "Organoid - Colon","Organoid - Colon", 
                                        "Organoid - Small Intestine","Organoid - Colon","Colon","Small Intestine", "Colon","Small Intestine",
                                        "Colon",  "Small Intestine","T-cell","Small Intestine" )


ggplot(Loadings_GEO_less, aes(PC1, -PC2, fill=tissue_plt, color=condition, shape=age))+geom_point(size=3)+theme_bw()+
  scale_fill_manual(values=c("#60B95D","#e9464b","#2b5826", "grey","#08306b", 
                             "cornflowerblue", "#871013"), name="Tissue")+
  scale_color_manual(values=c("black", "white"), name="Age Group")+th+
  scale_shape_manual(values = c(21, 22, 21))+
  guides(fill = guide_legend(override.aes=list(shape=21)),
         color = guide_legend(override.aes=list(shape=21)))

ggsave(here("DNAm/figs","PC1_PC2_GEO_organoid_fetalonly.pdf"), width = 8.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_GEO_organoid_fetalonly.jpeg"), width = 8.5, height = 6)




#' ## PCA with just organoid and blood
                # load("../../data/GEO_ibd_beta_multi_fetal_organoid.RData")
                # load("../../data/PCA_with_organoids.RData")

Loadings_GEO$condition<-Loadings_GEO$tissue_plt
Loadings_GEO$condition[grep("organoid",Loadings_GEO$tissue_plt)]<-"Organoid"
Loadings_GEO$condition[which(Loadings_GEO$condition!="Organoid")]<-"Primary"

keep_ids<-Loadings_GEO$array.id[grep("T-cell|small|organoid|large|SC|Blood|purified|organism|TI|AC",Loadings_GEO$tissue_plt)]
keep_ids<-keep_ids[which(!(keep_ids%in%Loadings_GEO$array.id[grep("Fetal|fetal|foetal|stomach",Loadings_GEO$age)]))]
keep_ids<-keep_ids[which(!(keep_ids%in%Loadings_GEO$array.id[grep("stomach",Loadings_GEO$tissue_plt)]))]
keep_ids<-gsub("X","",keep_ids)

GEO_ibd_beta_less<-GEO_ibd_beta_multi_fetal_organoid[,which(colnames(GEO_ibd_beta_multi_fetal_organoid)%in%c(keep_ids))]

pca_res_GEO_less <- prcomp(t(GEO_ibd_beta_less[complete.cases(GEO_ibd_beta_less),]))
Loadings_GEO_less<-as.data.frame(pca_res_GEO_less$x)
Loadings_GEO_less$array.id<-rownames(Loadings_GEO_less)


Loadings_GEO_meta<-Loadings_GEO[,452:456]

                  # save(Loadings_GEO_less,Loadings_GEO_meta , file="../../data/PCA_with_organoids_nofetalorextra.RData")
                  #
                  #
                  # load("../../data/PCA_with_organoids_nofetalorextra.RData")
Loadings_GEO_meta$array.id<-gsub("X","",Loadings_GEO_meta$array.id)
Loadings_GEO_less<-merge(Loadings_GEO_less, Loadings_GEO_meta, by="array.id")
unique(Loadings_GEO_less$tissue_plt)[order(unique(Loadings_GEO_less$tissue_plt))]

Loadings_GEO_less$condition<-Loadings_GEO_less$tissue_plt
Loadings_GEO_less$condition[grep("organoid",Loadings_GEO_less$tissue_plt)]<-"Organoid"
Loadings_GEO_less$condition[which(Loadings_GEO_less$condition!="Organoid")]<-"Primary"

Loadings_GEO_less$tissue_plt<-gsub("organism part |purified cell ","",Loadings_GEO_less$tissue_plt)
Loadings_GEO_less$tissue_plt<-as.factor(Loadings_GEO_less$tissue_plt)
levels(Loadings_GEO_less$tissue_plt)<-c("AC","AC","Whole Blood", "Large Intestine",
                                        "Organoid - Gastric Heterotopia", "Organoid - Rectum","Organoid - SC", "Organoid - TI",
                                        "SC", "SC", "Small Intestine","T-cell","TI" ,"TI" )

ggplot(Loadings_GEO_less, aes(PC1, PC2, fill=tissue_plt, color=study, size=study))+geom_point(shape=21)+theme_bw()+
  scale_fill_manual(values=c("#1a9850","#cb181d", "#60B95D","#004529",
                             "#004529","#004529",
                             "#08306b","#a6d96a","cornflowerblue", "#cb181d", "cornflowerblue"), name="Tissue")+
  th+scale_color_manual(values=c("black","black", "white"))+scale_size_manual(values=c(3,3,2))

ggsave(here("DNAm/figs","PC1_PC2_GEO_organoid.pdf"), width = 8.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_GEO_organoid.jpeg"), width = 8.5, height = 6)



#' ## highlight rescopes
Loadings_GEO_less$rescopeID<-sapply(1:nrow(Loadings_GEO_less), function(x){
  if(Loadings_GEO_less$array.id[x]%in%sampleinfo$array.id){
    sampleinfo$sample_ID[which(sampleinfo$array.id==Loadings_GEO_less$array.id[x])]
  }else{Loadings_GEO_less$array.id[x]}
})


ggplot(Loadings_GEO_less, aes(PC1, PC2, fill=tissue_plt, color=study))+geom_point(shape=21, aes(size=study))+theme_bw()+
  geom_line(aes(group=rescopeID))+
  scale_fill_manual(values=c("#1a9850","#cb181d", "#60B95D","#004529",
                             "#004529","#004529",
                             "#08306b","#a6d96a","cornflowerblue", "#cb181d", "cornflowerblue"), name="Tissue")+
  th+scale_color_manual(values=c("black","black", "white"))+scale_size_manual(values=c(3,3,2))

ggsave(here("DNAm/figs","PC1_PC2_GEO_organoid_rescopes.pdf"), width = 8.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_GEO_organoid_rescopes.jpeg"), width = 8.5, height = 6)
















#' ## More T cells
terms<-"T cell|t cell|CD4+"
tcell_meta<-meta[unique(c(grep(terms,meta$source_name_ch1),
                          grep(terms,meta$title),
                          grep("terms",meta$characteristics_ch1))) ,]
print(paste("Gut Samples: ",nrow(tcell_meta), "   Gut Studies: ", length(unique(tcell_meta$series_id)), sep=""))

tcell_meta<-tcell_meta[-grep("Guthrie|FFPE|adenocarcinoma|after treatment|tumor|thymus|Colon|cancer|oma|stem cells|HIV|cord blood|HCT116|Relapsing|Secondary|Old", tcell_meta$source_name_ch1),]
tcell_meta<-tcell_meta[-grep("cancer|cord blood|HIV|Lupus|lupus|12m|neonatal|Multiple sclerosis|CMV|fibroblast|GSE71955", tcell_meta$characteristics_ch1),]
#Disease studies not controls
tcell_meta<-tcell_meta[-grep("GSE89251|GSE56047|GSE59250|GSE71955", tcell_meta$series_id),]
#very few samples
fewsample<-names(tapply(tcell_meta$gsm,tcell_meta$series_id, function(x) length(unique(x)))[which(tapply(tcell_meta$gsm,tcell_meta$series_id, function(x) length(unique(x)))<10)])
tcell_meta<-tcell_meta[which(!(tcell_meta$series_id%in%fewsample)),]

print(paste("T Cell Samples: ",nrow(tcell_meta), "   T Cell Studies: ", length(unique(tcell_meta$series_id)), sep=""))





#' ### Going to collect GSE87640 GSE79237

GSE87640<-meta[grep("GSE87640", meta$series_id),]
#just t cells
GSE87640<-GSE87640[grep("CD4|CD8",GSE87640$source_name_ch1),]
GSE87640_beta<- as.data.frame(exprs(getGEO("GSE87640")[[1]]))
GSE87640_beta<-GSE87640_beta[which(colnames(GSE87640_beta)%in%GSE87640$gsm)] # 115 samples of CD4 or CD8 T cells from 59 people
dim(GSE87640_beta)
save(GSE87640_beta, file=paste(here("DNAm/data/public/"),"GSE87640_beta.RData", sep=""))


GSE79237<-meta[grep("GSE79237", meta$series_id),]
GSE79237_beta<- as.data.frame(exprs(getGEO("GSE79237")[[1]]))
GSE79237_beta<-GSE79237_beta[,which(colnames(GSE79237_beta)%in%GSE79237$gsm)] # 66 samples of CD4  T cells
dim(GSE79237_beta)
save(GSE79237_beta, file=paste(here("DNAm/data/public/"),"GSE79237_beta.RData", sep=""))


GSE50222_beta<-GSE50222_beta[which(rownames(GSE50222_beta)%in%rownames(ibd_combo)),]
ibd_combo_joint<-ibd_combo[which(rownames(ibd_combo)%in%rownames(GSE50222_beta)),]
GSE50222_beta<-GSE50222_beta[match(rownames(ibd_combo_joint),rownames(GSE50222_beta)),]
identical(rownames(ibd_combo_joint), rownames(GSE50222_beta))

GSE87640_beta<-GSE87640_beta[which(rownames(GSE87640_beta)%in%rownames(ibd_combo)),]
ibd_combo_joint<-ibd_combo[which(rownames(ibd_combo)%in%rownames(GSE87640_beta)),]
GSE87640_beta<-GSE87640_beta[match(rownames(ibd_combo_joint),rownames(GSE87640_beta)),]
identical(rownames(ibd_combo_joint), rownames(GSE87640_beta))
identical(rownames(ibd_combo_joint), rownames(ibd_combo_joint))

GSE79237_beta<-GSE79237_beta[which(rownames(GSE79237_beta)%in%rownames(ibd_combo)),]
ibd_combo_joint<-ibd_combo[which(rownames(ibd_combo)%in%rownames(GSE79237_beta)),]
GSE79237_beta<-GSE79237_beta[match(rownames(ibd_combo_joint),rownames(GSE79237_beta)),]
identical(rownames(ibd_combo_joint), rownames(GSE79237_beta))

identical(rownames(GSE50222_beta), rownames(GSE79237_beta))
identical(rownames(GSE50222_beta), rownames(GSE87640_beta))


#' ### Combat T cells
tcell_beta<-cbind(GSE50222_beta, GSE87640_beta, GSE79237_beta)
tcell_sampleinfo<-data.frame(sample_ID=c(colnames(GSE50222_beta),colnames(GSE87640_beta),colnames(GSE79237_beta)),
                             Tissue=c(rep("CD4+ T Cell", ncol(GSE50222_beta)),
                                      GSE87640$source_name_ch1, GSE79237$source_name_ch1),
                             Condition=c(rep("Primary",ncol(GSE50222_beta)),rep("Primary",ncol(GSE87640_beta)),rep("Primary",ncol(GSE79237_beta))),
                             study=c(rep("GSE50222",ncol(GSE50222_beta)),rep("GSE87640",ncol(GSE87640_beta)),rep("GSE79237",ncol(GSE79237_beta))))

identical(tcell_sampleinfo$sample_ID, colnames(tcell_beta))


tcell_sampleinfo$study<-as.factor(tcell_sampleinfo$study)

# # impute 0 and 1
tcell_beta[tcell_beta==0]<-0.01
tcell_beta[tcell_beta==1]<-0.99

# impute NA
imputeMedianv3<-function(x) apply(x, 1, function(x){x[is.na(x)]<-median(x, na.rm=T); x}) #impute with row mean
tcell_beta<-t(imputeMedianv3(tcell_beta))

Mval<-function(beta) log2(beta/(1-beta))
edata = apply(tcell_beta, 1, Mval) # need mvalues for combat
edata = as.data.frame(edata)
edata = t(edata)



#mod = model.matrix(~as.factor(Tissue), data=tcell_sampleinfo) # can not protect cause confounded
batch = tcell_sampleinfo$study
combat_tcell_mval = ComBat(dat=edata, batch=batch, mod=NULL, par.prior=TRUE)

#Back to betas
betas<-function(M) 2^M/((2^M)+1)
combat_tcell_Beta = apply(combat_tcell_mval, 1, betas) # need mvalues for combat
combat_tcell_Beta = as.data.frame(combat_tcell_Beta)
combat_tcell_Beta = t(combat_tcell_Beta)
combat_tcell_Beta<-as.data.frame(combat_tcell_Beta)

save(combat_tcell_Beta, tcell_sampleinfo, file=paste(here("DNAm/data/public/"),"tcells_combatted_GEO.Rdata", sep=""))


## bring in organoids
sampleinfo_organoid <- read.table(here("DNAm/data/public/organoid/E-MTAB-4957.sdrf.txt"), header=T, sep="\t")
sampleinfo_organoid<-sampleinfo_organoid[,c(1:15,30:34,38:41)]
sampleinfo_organoid$sentrix<-sapply(1:nrow(sampleinfo_organoid), function(x) strsplit(as.character(sampleinfo_organoid$Assay.Name)[x], "_")[[1]][1])

sampleinfo_organoid<-sampleinfo_organoid[-grep("Fetal|fetal|foetal|stomach",sampleinfo_organoid$Factor.Value.developmental.stage.),]
organoid_beta<-organoid_beta[,colnames(organoid_beta)%in%sampleinfo_organoid$Assay.Name]
identical(colnames(organoid_beta), as.character(sampleinfo_organoid$Assay.Name))

identical(rownames(organoid_beta), rownames(combat_tcell_Beta))
identical(rownames(organoid_beta), rownames(ibd_combo))

GEO_ibd_beta<-cbind(ibd_combo, organoid_beta,combat_tcell_Beta)

#'## PCA with just organoid and T cell
pca_res_GEO <- prcomp(t(GEO_ibd_beta[complete.cases(GEO_ibd_beta),]))
vars <- pca_res_GEO$sdev^2
Importance<-vars/sum(vars)

print("the variance of PC1 and PC2 are", round(Importance[1],digits =2), round(Importance[2],digits =2))

Loadings_GEO<-as.data.frame(pca_res_GEO$x)
Loadings_GEO$array.id<-rownames(Loadings_GEO)

GSE87640<-GSE87640[which(GSE87640$gsm%in%colnames(GSE87640_beta)),]
GSE79237<-GSE79237[which(GSE79237$gsm%in%colnames(GSE79237_beta)),]

sampleinfo_organoid$tissue<-paste(sampleinfo_organoid$Characteristics.biosource.type., sampleinfo_organoid$Factor.Value.sampling.site.)

Loadings_GEO_meta<-data.frame(sample_ID=c(colnames(ibd_combo),colnames(organoid_beta),colnames(GSE50222_beta),colnames(GSE87640_beta),colnames(GSE79237_beta)),
                              Tissue=c(as.character(sampleinfo$sample.site_grouped), as.character(sampleinfo_organoid$Factor.Value.sampling.site.),rep("T Cell", ncol(GSE50222_beta)),
                                       GSE87640$source_name_ch1, GSE79237$source_name_ch1),
                              Condition=c(rep("Primary",ncol(ibd_combo)),as.character(sampleinfo_organoid$Characteristics.biosource.type.),rep("Primary",ncol(GSE50222_beta)),rep("Primary",ncol(GSE87640_beta)),rep("Primary",ncol(GSE79237_beta))),
                              study=c(rep("IBD",ncol(ibd_combo)),rep("E-MTAB-4957",ncol(organoid_beta)),rep("GSE50222",ncol(GSE50222_beta)),rep("GSE87640",ncol(GSE87640_beta)),rep("GSE79237",ncol(GSE79237_beta))))
Loadings_GEO_meta$Tissue<-as.factor(Loadings_GEO_meta$Tissue)
levels(Loadings_GEO_meta$Tissue)<-c("Colon","CD4+ T Cell","CD8+ T Cell","Colon","Other GI","Small Intestine","CD4+ T Cell", "Colon","Colon","Other GI","CD4+ T Cell","Small Intestine")
Loadings_GEO_meta$Condition<-as.factor(Loadings_GEO_meta$Condition)
levels(Loadings_GEO_meta$Condition)<-c("Primary","Organoid","Primary","Primary")
            #
            # save(Loadings_GEO,Loadings_GEO_meta , file="../../data/PCA_with_organoids_extraTCells.RData")
            #
            #
            # load("../../data/PCA_with_organoids_extraTCells.RData")
Loadings_GEO$array.id<-gsub("X","",rownames(Loadings_GEO))
Loadings_GEO<-merge(Loadings_GEO, Loadings_GEO_meta, by.x="array.id", by.y="sample_ID")

Loadings_GEO$label<-sapply(1:nrow(Loadings_GEO), function(x) {
  if(Loadings_GEO$Condition[x]=="Organoid"){paste(Loadings_GEO$Condition[x]," - ", Loadings_GEO$Tissue[x])}else{as.character(Loadings_GEO$Tissue[x])}})
Loadings_GEO$label<-as.factor(Loadings_GEO$label)


ggplot(Loadings_GEO, aes(PC1, PC2, fill=label, color=study))+geom_point(shape=21, size=3)+theme_bw()+
  scale_fill_manual(values=c("#871013","#e9464b","#60B95D","#2b5826", "grey", "#08306b","cornflowerblue"), name="Tissue")+
  th+scale_color_manual(values=c("white","white", "white", "white", "black")) + xlab("PC1 (71%)") + ylab("PC2 (6%)")


ggsave(here("DNAm/figs","PC1_PC2_GEO_organoid_moreTcells.pdf"), width = 8.5, height = 6)
ggsave(here("DNAm/figs/jpeg","PC1_PC2_GEO_organoid_moreTcells.jpeg"), width = 8.5, height = 6)






#'## R Session Info
sessionInfo()
