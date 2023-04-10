#'---
#'title: Annotation of EPIC CpGs to genes
#'author: Rachel Edgar
#'date: "`r Sys.Date()`"
#'---
#'

suppressMessages(library(here))

options(stringsAsFactors = FALSE)
options(scipen = 999)


#' system arguments
args <- commandArgs(trailingOnly = TRUE)
cpg_1<-as.numeric(args[1])
if(cpg_1>860000){cpg_2<-865859}else{cpg_2<-cpg_1+9999}



load(here("DNAm/data","ensembl_transcripts.RData"))

#' CpG Coordiantes from illumina annotation
anno_EPIC<-read.csv(here("data", "MethylationEPIC_v-1-0_B4.csv"), skip=7)
anno_EPIC_minimal<-anno_EPIC[,c("IlmnID","CHR", "MAPINFO")]
anno_EPIC_minimal<-anno_EPIC_minimal[which(!(is.na(anno_EPIC_minimal$MAPINFO))),]


CpG_to_gene<-lapply(cpg_1:cpg_2, function(x){

  chr_transcripts<- ensembl_transcripts[which(ensembl_transcripts$Chromosome.scaffold.name==anno_EPIC_minimal$CHR[x]),]
  
  genes<-lapply(1:nrow(chr_transcripts), function(y){
    if(anno_EPIC_minimal$MAPINFO[x] >= chr_transcripts$total_start[y] & anno_EPIC_minimal$MAPINFO[x] <= chr_transcripts$total_end[y]){
      if(anno_EPIC_minimal$MAPINFO[x] >= chr_transcripts$promoter_start[y] & anno_EPIC_minimal$MAPINFO[x] <= chr_transcripts$promoter_end[y]){feature<-"promoter"}else{
        if(anno_EPIC_minimal$MAPINFO[x] >= chr_transcripts$intragenic_start[y] & anno_EPIC_minimal$MAPINFO[x] <= chr_transcripts$intragenic_end[y]){feature<-"intragenic"}else{
          if(anno_EPIC_minimal$MAPINFO[x] >= chr_transcripts$UTR3_start[y] & anno_EPIC_minimal$MAPINFO[x] <= chr_transcripts$UTR3_end[y]){feature<-"3'UTR"}}
      }
      cbind(anno_EPIC_minimal[x,], chr_transcripts[y,c(1,4,9,10)], data.frame(CpG_in=feature))}else{
      }
  })
  
  genes<-do.call(rbind, genes)})

CpG_to_gene<-do.call(rbind, CpG_to_gene)


#' keep features that are different between transcripts but ditch the transcript IDs
CpG_to_gene_less<-unique(CpG_to_gene[,-5])

write.csv(CpG_to_gene_less, file=paste(here("DNAm/data/annotation_split"),"/EPIC_annotation_gene_",cpg_1,".csv", sep=""))
